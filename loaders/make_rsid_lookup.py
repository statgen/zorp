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
import argparse
import gzip
import itertools
import os
import re
import struct
import sys
import typing as ty

from fastnumbers import int
import lmdb
import msgpack


RSID_CAPTURE = re.compile(r'RS=(\d+);?')

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


def make_file_iterator(handle) -> ty.Iterator[ty.Tuple[str, int, str, str, int]]:
    """Parse the set of lines for an open file handle (text, gz, etc)"""
    for row in handle:
        if row.startswith('#'):
            # Skip header rows
            continue

        if not row.startswith('NC_'):
            # After that, only parse the variants labeled "NC", not the "unplaced scaffold" items
            raise StopIteration

        yield line_parser(row)


def make_group_iterator(file_iterator) -> ty.Iterator[ty.Tuple[str, int, dict]]:
    """Create an iterator that returns all possible ref/alt : rsid pairs for that position"""

    chrom = None
    for position, position_records in itertools.groupby(file_iterator, key=lambda x: x[1]):
        # assumes that position is the only thing that will change, else use chr/pos: (x[0], x[1])):
        position_contents = {}
        # Then push final lookup for that chr:pos into the database
        for record in position_records:
            # The same line can indicate one or more ref/alts
            chrom, pos, ref, alt, rsid = record
            ref_alts = [f'{ref_option}/{alt_option}'
                        for alt_option in alt.split(',')
                        for ref_option in ref.split('.')]

            for key in ref_alts:
                position_contents[key] = rsid

        yield chrom, position, position_contents


def parse_args():
    parser = argparse.ArgumentParser(description="Create a fast lookup for rsIDs")
    parser.add_argument('source',
                        help='The dbSNP VCF (in gzip format) that will be used to generate the lookup')
    parser.add_argument('out', help='Output filename for the resulting database')
    parser.add_argument('--force', action='store_true',
                        help='If the file exists, force it to be overwritten')
    parser.add_argument('--n_chroms', type=int, default=25,
                        help="Number of chromosomes in the input file")
    return parser.parse_args()


def make_databases(env) -> dict:
    known = [
        '1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
        '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
        '21', '22', 'X', 'Y', 'MT'
    ]
    return {k: env.open_db(bytes(k, "utf8"), integerkey=True) for k in known}


def main(source_fn: str, out_fn: str, force=False, n_chroms=25):
    if os.path.exists(out_fn):
        # LMDB would happily append / replace keys, but we want to create these files once and keep file size small
        print('The requested output file already exists.')
        if force:
            print('It will be deleted and recreated')
            os.remove(out_fn)
        else:
            sys.exit(1)

    # To reduce disk usage, each chromosome is stored internally as a separate database (the same chrom key info isn't
    #   stored redundantly)
    env = lmdb.open(out_fn, subdir=False, max_dbs=n_chroms, map_size=25e9)
    db_handles = make_databases(env)

    with env.begin(write=True) as txn, gzip.open(source_fn, "rt") as vcf:
        all_lines = make_file_iterator(vcf)
        group_iterator = make_group_iterator(all_lines)

        for chrom, position, position_contents in group_iterator:
            # Value is not a primitive; serialize it efficiently
            key = struct.pack('I', position)
            value = msgpack.packb(position_contents, use_bin_type=True)
            txn.put(key, value, db=db_handles[chrom])


if __name__ == "__main__":
    args = parse_args()
    main(args.source, args.out,
         force=args.force, n_chroms=args.n_chroms)
