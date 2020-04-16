"""
Convert SNP information to gene, rsid, or other useful annotations
"""
import struct
import typing as ty

import lmdb
import msgpack

from . import assets, exceptions


class SnpToRsid:
    """Convert SNP coordinates to RSID information"""
    def __init__(self, path: str, *, num_chroms: int = 25, test=False):
        if not path:
            raise exceptions.ConfigurationException('Must provide a path to the lookup file')

        if path == 'GRCh37' or path == 'GRCh38':
            record_type = 'snp_to_rsid'
            if test:
                record_type += '_test'
            path = assets.manager.locate(record_type, genome_build=path)

        self.env = lmdb.open(path, subdir=False, max_dbs=num_chroms, readonly=True)
        self.db_handles = {}  # type: dict

    def __call__(self, chrom: str, pos: int, ref: str, alt: str) -> ty.Union[int, None]:
        """
        Look up the specified SNP in the database, and handle any translation of how results are stored in this
            specific file format.
        """
        if chrom not in self.db_handles:
            self.db_handles[chrom] = self.env.open_db(bytes(chrom, 'utf8'), integerkey=True)

        db = self.db_handles[chrom]

        key = struct.pack('I', int(pos))
        with self.env.begin(buffers=True) as txn:
            res = txn.get(key, db=db)
            if res:
                res = msgpack.unpackb(res, use_list=False)
                res = res.get('{}/{}'.format(ref, alt))
        return res

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
