#!/usr/bin/env python3
"""
Generate an LMDB database optimized for chr:pos:ref:alt --> rsid lookups

Sample usage:
    ./make_rsid_lookups.py dbSNPvcf.gz dbSNPb153.lmbd --force


Each database file is structured as one file, multiple databases (one DB per chromosome)
Each per-chromosome DB has one key per (int)position, and each value is a hash of possible {"ref/alt":rsid} matches.
Positions and rsids are represented as integers to save space in an effort to stay below virtual memory limits.

The final database format (one database per chromosome) is:

{  pos:  { 'ref/alt1': 1 , 'ref/alt2': 2 }  (where rs1 and rs2 are represented as integers, "1" and "2")
"""
import gzip
import itertools
import os
import re
import struct
import tempfile
import typing as ty
import urllib.request

from fastnumbers import int
from filefetcher.manager import BuildTask
import lmdb
import msgpack
import pysam


RSID_CAPTURE = re.compile(r'RS=(\d+);?')

# The new dbSNP format uses refseq identifiers. Building a lookup based on human friendly chrom/pos/ref/alt
#  requires converting human identifiers to and from refseq names
VERSIONLESS_CHROMS = {
    'NC_000001': '1',
    'NC_000002': '2',
    'NC_000003': '3',
    'NC_000004': '4',
    'NC_000005': '5',
    'NC_000006': '6',
    'NC_000007': '7',
    'NC_000008': '8',
    'NC_000009': '9',
    'NC_000010': '10',
    'NC_000011': '11',
    'NC_000012': '12',
    'NC_000013': '13',
    'NC_000014': '14',
    'NC_000015': '15',
    'NC_000016': '16',
    'NC_000017': '17',
    'NC_000018': '18',
    'NC_000019': '19',
    'NC_000020': '20',
    'NC_000021': '21',
    'NC_000022': '22',
    'NC_000023': 'X',
    'NC_000024': 'Y',
    'NC_012920': 'MT',
}


def make_chrom_to_contigs(tabix_file: pysam.TabixFile) -> dict:
    """
    In order to address a tabix file by position, we need the exact contig name. dbSNP builds use a versioned
        identifier that changes across builds and versions.

    This generates a lookup for converting human friendly chromosome names to things that can be used by a
        particular dbSNP file. It automatically adjusts for small build differences and excludes non-chrom contigs
    """
    return {
        # Eg '1': 'NC_000001.10',
        VERSIONLESS_CHROMS[c.split('.', 1)[0]]: c
        for c in tabix_file.contigs if c.startswith('NC')
    }


def fetch_regions_sequentially(source_fn: str, sample_regions: ty.Tuple[str, int, int]) -> ty.Iterator[str]:
    """Instead of loading an entire gzip file, create an iterator over a set of regions"""
    source = pysam.TabixFile(source_fn)
    human_chrom_to_tabix = make_chrom_to_contigs(source)

    for region in sample_regions:
        chrom, start, end = region
        dbsnp_chrom = human_chrom_to_tabix[chrom]

        for line in source.fetch(dbsnp_chrom, start, end):
            yield line


def line_parser(row) -> ty.Tuple[str, int, str, str, int]:
    """For new dbSNP format, builds 152+"""
    fields = row.split()
    # the new dbSNP format uses refseq ids + version; convert these to human-readable chromosome names
    chrom = VERSIONLESS_CHROMS[fields[0].split('.')[0]]
    pos = int(fields[1])

    ref = fields[3]
    alt = fields[4]

    # Get the RSID from the VCF info field, in case the id column is ambiguous for some reason
    rsid = int(RSID_CAPTURE.search(fields[7]).group(1))

    return (chrom, pos, ref, alt, rsid)


def make_file_iterator(handle: ty.Iterable) -> ty.Iterator[ty.Tuple[str, int, str, str, int]]:
    """
    Parse the set of lines for some iterator (eg over contents of an open text, gz, etc file), and only give back the
        ones relevant to our chrom/pos/ref/alt use case
    """
    for row in handle:
        if row.startswith('#'):
            # Skip header rows
            continue

        if not row.startswith('NC_'):
            # After that, only parse the variants labeled "NC", not the "unplaced scaffold" items
            raise StopIteration

        yield line_parser(row)


def make_group_iterator(file_iterator) -> ty.Iterable[ty.Tuple[str, int, dict]]:
    """Create an iterator that returns all possible ref/alt : rsid pairs for that position"""

    chrom = None
    for position, position_records in itertools.groupby(file_iterator, key=lambda x: x[1]):
        # assumes that position is the only thing that will change, else use chr/pos: (x[0], x[1])):
        position_contents = {}
        # Then push final lookup for that chr:pos into the database
        for record in position_records:
            # The same line can indicate one or more ref/alts
            chrom, pos, ref, alt, rsid = record
            ref_alts = ['{}/{}'.format(ref_option, alt_option)
                        for alt_option in alt.split(',')
                        for ref_option in ref.split('.')]

            for key in ref_alts:
                position_contents[key] = rsid

        yield chrom, position, position_contents


def make_databases(env) -> dict:
    known = [
        '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
        '21', '22', 'X', 'Y', 'MT'
    ]
    return {k: env.open_db(bytes(k, "utf8"), integerkey=True) for k in known}


def main(source_fn: str, out_fn: str, n_chroms=25, sample_regions=None):
    """Perform this task in isolation for any dbsnp file"""
    if os.path.exists(out_fn):
        # LMDB would happily append / replace keys, but we want to create these files once and keep file size small
        raise Exception('The requested output file already exists.')

    # To reduce disk usage, each chromosome is stored internally as a separate database (the same chrom key info isn't
    #   stored redundantly)
    env = lmdb.open(out_fn, subdir=False, max_dbs=n_chroms, map_size=25e9)
    db_handles = make_databases(env)

    if sample_regions:
        iterator = fetch_regions_sequentially(source_fn, sample_regions)
    else:
        iterator = gzip.open(source_fn, "rt")

    with env.begin(write=True) as txn:
        all_lines = make_file_iterator(iterator)
        group_iterator = make_group_iterator(all_lines)

        for chrom, position, position_contents in group_iterator:
            # Value is not a primitive; serialize it efficiently
            key = struct.pack('I', position)
            value = msgpack.packb(position_contents, use_bin_type=True)
            txn.put(key, value, db=db_handles[chrom])


class MakeSnpToRsid(BuildTask):
    """A packaged filefetcher build task that also downloads the necessary input files"""
    def __init__(self, genome_build, sample_regions=None):
        self.genome_build = genome_build
        self.regions = sample_regions or None

    def get_assets(self):
        # Download the appropriate dbsnp file for a given genome build if does not exist
        # Uses a system temp directory rather than the build folder in case file is needed between runs

        # FIXME: The directory structure of new dbSNP is still too new to know what will change between releases, so
        #  we don't yet have a generic way to get build information without downloading the files. For
        #  now this is a hardcoded reference.
        if self.genome_build == 'GRCh37':
            source_url = 'ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.25.gz'
        elif self.genome_build == 'GRCh38':
            source_url = 'ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.38.gz'
        else:
            raise Exception('Unknown build')

        dest_fn = os.path.join(tempfile.gettempdir(), os.path.basename(source_url))

        if not os.path.exists(dest_fn):
            # Download assets to a tempfile directory
            urllib.request.urlretrieve(
                source_url,
                dest_fn
            )
            # And also Tabix index
            urllib.request.urlretrieve(
                source_url + '.tbi',
                dest_fn + '.tbi'
            )
        return dest_fn

    def build(self, manager, item_type: str, build_folder: str, **kwargs):
        # FIXME: Hardcode dbsnp build until we build a more reliable way to get this info.
        dbsnp_build = 'b153'

        source_fn = self.get_assets()
        print('Building to: ', build_folder)
        dest_fn = os.path.join(
            build_folder,
            '{}_{}_{}.lmdb'.format(item_type, self.genome_build, dbsnp_build)
        )

        main(source_fn, dest_fn, sample_regions=self.regions)

        return dest_fn, {'dbsnp_build': dbsnp_build}
