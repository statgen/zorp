#!/usr/bin/env python
"""
Look up the RSID associated with a variant, based on a tabixed copy of dbSNP
"""
import typing as ty

import pysam

from zorp import lookups, parsers, sniffers

from loaders.make_rsid_lookup import make_group_iterator, make_file_iterator

# Our datasets use human-readable chromosomes, which must be converted to exactly the format used in the tabix
#   index provided by NCBI (the .x suffix is a version string which may vary across builds; this is from b153)
# TODO In a production script, we'd control for version numbers
CHROM_TO_TABIX = {
    '1': 'NC_000001.10',
    '10': 'NC_000010.10',
    '11': 'NC_000011.9',
    '12': 'NC_000012.11',
    '13': 'NC_000013.10',
    '14': 'NC_000014.8',
    '15': 'NC_000015.9',
    '16': 'NC_000016.9',
    '17': 'NC_000017.10',
    '18': 'NC_000018.9',
    '19': 'NC_000019.9',
    '2': 'NC_000002.11',
    '20': 'NC_000020.10',
    '21': 'NC_000021.8',
    '22': 'NC_000022.10',
    '3': 'NC_000003.11',
    '4': 'NC_000004.11',
    '5': 'NC_000005.9',
    '6': 'NC_000006.11',
    '7': 'NC_000007.13',
    '8': 'NC_000008.10',
    '9': 'NC_000009.11',
    'MT': 'NC_012920.1',
    'X': 'NC_000023.10',
    'Y': 'NC_000024.9'
}


class LookupRsidsTabix:
    """
    Find RSIDs that match a specified file. This is a tabix-based design and inherently assumes that each lookup
        performed will be sequential (after the previous one)
    """
    def __init__(self, source_fn: str, max_segment: int = 50000):
        # Define the data to iterate over
        self._tabix = pysam.TabixFile(source_fn)
        self._current_reader = None  # type: ty.Iterator[str]

        # Store the "current" row (sometimes our dbsnp reader will get AHEAD of the sumstats reader,
        #   or be asked for the same position twice in a row, eg multiallelic variants)
        self._current_row = None

        # Optimization: actual dbSNP is a very large file, and we want to only iterate over the parts near variants
        #   in the sumstats file
        # To avoid excessive disk seeks, we'll prefetch 50kb regions at a time

        # TODO: Investigate tabix file
        #   1. Get the chromosome names from the tabix index provided by NIH (remember, columns are a bit weird)
        #   2. Use chromosome ranges to grab tabix data (requires pysam instead of gzip reader)

        # Identify where the last tabix query was run
        self._last_query_chrom = None  # In human readable format so we can quickly compare to the current row
        self._last_query_pos = 0
        self._max_segment = max_segment  # type: int

    def make_reader(self, chrom, start_pos):
        """Make a tabix reader from the data for the specified region"""
        # TODO in a production tool we'd validate the chromosome to ensure it was known to dbsnp
        dnsnp_chrom = CHROM_TO_TABIX[chrom]
        segment = self._tabix.fetch(dnsnp_chrom, start_pos, start_pos + self._max_segment)
        all_lines = make_file_iterator(segment)

        # Update current reader position
        self._current_reader = iter(make_group_iterator(all_lines))
        self._current_row = next(self._current_reader)

        self._last_query_chrom = chrom
        self._last_query_pos = start_pos

    def _advance_current_reader(self, target_chrom, target_pos):
        # 1. Decide whether to advance the current reader (a 50kb chunk of data) or seek to a new section of the
        #   file on disk
        if self._last_query_chrom != target_chrom or target_pos > (self._last_query_pos + self._max_segment):
            self.make_reader(target_chrom, target_pos)

        # 2. Advance the reader to the desired position
        chrom, pos, _ = self._current_row

        if chrom == target_chrom and pos >= target_pos:
            # Don't advance the dbSNP reader any further- we are already at or past where we need to be
            return self._current_row

        # If the dbSNP reader is lagging the gwas reader, then we need to advance
        while chrom != target_chrom or pos < target_pos:
            try:
                self._current_row = next(self._current_reader)
            except StopIteration:
                # If we hit the end of the dbSNP file, exit the loop and never let the iterator advance again
                break

            chrom, pos, _ = self._current_row

    def __call__(self, chrom: str, pos: int, ref: str, alt: str) -> ty.Union[int, None]:
        # Look up the rsid associated with a current SNP.
        self._advance_current_reader(chrom, pos)

        snp_chrom, snp_pos, ref_alt_options = self._current_row
        if snp_chrom != chrom or snp_pos != pos:
            # Ensure that the dbSNP reader is at the correct position!
            return None
        return ref_alt_options.get('{}/{}'.format(ref, alt), None)

    def __del__(self):
        self._tabix.close()


if __name__ == '__main__':
    ####
    # Local development
    # gwas_reader = sniffers.guess_gwas_standard(
    #     '/Users/abought/code/personal/zorp/loaders/benchmarks/data/summary_stats_1k_lines.tab'
    # )
    # TODO: Make an LMDB version based on just part of the tabix file
    # Eg: sumstats = text or gzip file
    # rsid_finder = lookups.FindRsid('/home/abought/dbsnp/b153/dbSNP_grch37_b153.lmdb')

    # rsid_finder = LookupRsidsTabix(
    #     '/Users/abought/code/personal/zorp/loaders/benchmarks/data/dbsnp_b153_fragment_for_sumstats_benchmarks.gz'
    # )

    ###
    # Options useful on the cluster (full dataset)
    gwas_reader = sniffers.guess_gwas_standard('/home/abought/code/zorp/loaders/benchmarks/data/summary_stats.gz')
    rsid_finder = LookupRsidsTabix('/home/abought/dbsnp/b153/GCF_000001405.25.gz')
    # rsid_finder = lookups.FindRsid('/home/abought/dbsnp/b153/dbSNP_grch37_b153.lmdb')

    # Perform lookups and track results
    gwas_reader.add_lookup('rsid', lambda variant: variant.rsid or rsid_finder(variant.chrom, variant.pos, variant.ref,
                                                                               variant.alt))
    gwas_reader.write('deleteme.gz')
