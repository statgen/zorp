"""
Generate an LMDB database optimized for chr:pos:ref:alt --> rsid lookups

Each database file is structured as one file, multiple databases (one DB per chromosome)
Each per-chromosome DB has one key per position, and each value is a hash of possible {refalt:rsid} matches.
Chromosomes and rsids are represented as integers to save space in an effort to stay below virtual memory limits.

Known chromosomes (according to dbsnp): 25 in all. It seems we need to know this number in advance to use LMDB.
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, X, Y, MT]
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

    # To reduce space usage, each chromosome is stored internally as a separate database (the same chrom key info isn't
    #   stored redundantly)
    env = lmdb.open(out_fn, subdir=False, max_dbs=n_chroms, map_size=2e9)
    db_handles = make_databases(env)

    with env.begin(write=True) as txn, gzip.open(source_fn, "rt") as vcf:
        cur_pos = None
        cur_data = {}
        for line in vcf:
            chrom, pos, vid, ref, alt = line.split("\t")[0:5]
            rsid = vid.split(":")[0]

            if not rsid.startswith("rs"):
                continue

            # Some file size reduction by removing extra characters and encoding rsid as int
            rsid = int(rsid.replace('rs', ''))
            if cur_pos is None or pos == cur_pos:
                # A dbSNP file can have several lines for same position. Consolidate into a single key.
                cur_data[f'{ref}/{alt}'] = rsid
            else:
                # Push any existing lookup to the database and reset cur_data
                key = struct.pack('I', int(pos))
                value = msgpack.packb(cur_data, use_bin_type=True)  # Value is not a primitive; serialize it efficiently
                txn.put(key, value, db=db_handles[chrom])
                cur_data = {}

            cur_pos = pos

        # Make sure that the last line of the file gets added to the database (TODO: dry with logic above)
        # key = int(pos).to_bytes(4, 'big')  # 4 bytes = ~4 billion positions (unsigned int)
        value = msgpack.packb(cur_data, use_bin_type=True)  # Value is not a primitive; serialize it efficiently
        txn.put(key, value, db=db_handles[chrom], append=True)


if __name__ == "__main__":
    args = parse_args()
    main(args.source, args.out,
         force=args.force, n_chroms=args.n_chroms)
