"""
Convert SNP information to gene, rsid, or other useful annotations
"""
import struct

import lmdb
import msgpack


class FindRsid:
    """Convert SNP coordinates to RSID information"""
    def __init__(self, path: str, *, num_chroms: int = 25):
        if path == 'GRCh37' or path == 'GRCh38':
            raise NotImplementedError('Pre-defined builds and data cache directory not yet implemented')

        # TODO: Hardcoded, machine-specific default dir. Replace!
        path = path or '/Users/abought/code/personal/zorp/tests/data/rsid-build37-segment.lmdb'

        self.env = lmdb.open(path, subdir=False, max_dbs=num_chroms, readonly=True)
        self.db_handles = {}

    def __call__(self, chrom: str, pos: int, ref: str, alt: str):
        """
        Look up the specified SNP in the database, and handle any translation of how results are stored in this
            specific file format.
        """
        if chrom not in self.db_handles:
            self.db_handles[chrom] = self.env.open_db(bytes(chrom, 'utf8'), integerkey=True)

        db = self.db_handles[chrom]

        pos = struct.pack('I', int(pos))
        with self.env.begin(buffers=True) as txn:
            res = txn.get(pos, db=db)
            if res:
                res = msgpack.unpackb(res, encoding="utf8",use_list=False)
                res = res.get(f'{ref}/{alt}')
        return f'rs{res}' if res else None

    def known_chroms(self):
        with self.env.begin() as txn:
            # In this particular data structure, all child databases are listed as keys in the default DB
            cursor = txn.cursor()
            return [k for k, v in cursor]

    def __del__(self):
        # Ensure the database environment is closed when the reader is no longer needed
        try:
            self.env.close()
        except Exception:
            pass
