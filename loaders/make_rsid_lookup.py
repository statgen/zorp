"""
Generate an LMDB database optimized for chr:pos:ref:alt --> rsid lookups

Each database file is structured as one file, multiple databases (one DB per chromosome)
Each per-chromosome DB has one key per position, and each value is a hash of possible {"ref/alt":rsid} matches.
Positions and rsids are represented as integers to save space in an effort to stay below virtual memory limits.

Known chromosomes (according to dbsnp): 25 in all. It seems we need to know this number in advance to use LMDB.
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT]


Useful metrics:
- 4 byte unsigned int : max ~4bln. Even with TOPMED growth, probably won't hit that for a while.
- Longest human chromosome ~250Mb (chromosome 1), fits in 4 byte unsigned int
- Number of known rsids, Oct 2019: 660,146,231 (each fits in 4 byte unsigned int) (from wc -l)


Expected number entries per chromosome (!= #SNPS because some SNPs are multiallelic):
    zcat human9606_vcf_00-All.vcf.gz | tail -n +58 | cut -f1 | uniq -c

51651989 1
55499999 2
45421886 3
43657384 4
40973683 5
38315678 6
36567772 7
34770017 8
28774299 9
30525060 10
31275303 11
30324935 12
22388720 13
20383978 14
19100191 15
21009708 16
18595322 17
17708039 18
14227302 19
14550131 20
8717639 21
9064959 22
26166816 X
 473067 Y
   2297 MT


Consolidated LMDB: 634 M entries vs 660 M rsid lines (~3.8% multiallelic snps. A 2015 paper says "2.3%" was an old and
 likely low estimate, so this is not unreasonable.
"""
import argparse
import gzip
import os
import struct
import sys

import lmdb
import msgpack


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


def make_databases(env):
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
        prev_chrom = None
        prev_pos = None
        cur_data = {}
        for line in vcf:
            if line.startswith('#'):
                # Skip header rows
                continue

            chrom, pos, vid, ref, alt = line.split("\t")[0:5]
            rsid = vid.split(":")[0]

            if not rsid.startswith("rs"):
                continue

            # A dbSNP file can have several lines for same position. For each line, decide whether to BUILD a lookup
            # (same as previous line), or WRITE it to disk (a new position)

            if (prev_pos is not None and prev_pos != pos) or (prev_chrom is not None and prev_chrom != chrom):
                # This is a new position. Write the data from the old position before continuing
                # Push any existing lookup to the database and reset cur_data
                key = struct.pack('I', int(prev_pos))
                value = msgpack.packb(cur_data, use_bin_type=True)  # Value is not a primitive; serialize it efficiently
                txn.put(key, value, db=db_handles[prev_chrom])
                cur_data = {}

            # Some file size reduction by removing extra characters and encoding rsid as int
            rsid = int(rsid.replace('rs', ''))
            cur_data['{}/{}'.format(ref, alt)] = rsid

            prev_pos = pos
            prev_chrom = chrom

        # After reading all lines, ensure that the last item of data is also written to output
        key = struct.pack('I', int(prev_pos))
        value = msgpack.packb(cur_data, use_bin_type=True)  # Value is not a primitive; serialize it efficiently
        txn.put(key, value, db=db_handles[prev_chrom])


if __name__ == "__main__":
    args = parse_args()
    main(args.source, args.out,
         force=args.force, n_chroms=args.n_chroms)
