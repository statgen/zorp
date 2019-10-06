"""
Convert SNP information to gene, rsid, or other useful annotations
"""
import sys

import lmdb
import msgpack


class FindRsid:
    """Convert SNP coordinates to RSID information"""
    def __init__(self, path: str, *, num_chroms: int= 25):
        if path == 'GRCh37' or path == 'GRCh38':
            raise NotImplementedError('Pre-defined builds and data cache directory not yet implemented')

        # TODO: Hardcoded, machine-specific default dir. Replace!
        path = path or '/Users/abought/code/personal/zorp/tests/data/rsid-build37-segment.lmdb'

        self.env = lmdb.open(path, subdir=False, max_dbs=num_chroms)

    def __call__(self, chrom: str, pos: int, ref: str, alt: str):
        """
        Look up the specified SNP in the database, and handle any translation of how results are stored in this
            specific file format.
        """
        db = self.env.open_db(bytes(chrom, 'utf8'))  # always retrieves handle to same db
        pos = pos.to_bytes(4, 'little')  # FIXME: pos not found in test db
        with self.env.begin(buffers=True) as txn:
            res = txn.get(pos, db=db)
            if res:
                res = msgpack.unpackb(res,encoding="utf8",use_list=False)
                res = res.get(f'{ref}/{alt}')
        return f'rs{res}' if res else None

    def __del__(self):
        # Ensure the database environment is closed when the reader is no longer needed
        try:
            self.env.close()
        except Exception:
            pass
