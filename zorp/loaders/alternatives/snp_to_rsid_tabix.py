#!/usr/bin/env python
"""
Look up the RSID associated with a variant, based on a tabixed copy of dbSNP

This script is provided for internal development reference. It is not a fully-fledged lookup utility and should not
    be relied on for production use. Some filenames/ chrom versions may be hard-coded, etc.

TODO: Small changes to data pipeline have been made; verify that the `LookupRsidsTabix` still works
"""
import typing as ty

import pysam

from zorp import sniffers

from zorp import lookups
from zorp.loaders.make_rsid_lookup import make_group_iterator, make_file_iterator, make_chrom_to_contigs


class LookupRsidsTabix:
    """
    Find RSIDs that match a specified file. This is a tabix-based design and inherently assumes that each lookup
        performed will be sequential (after the previous one)
    """
    def __init__(self, source_fn: str, max_segment: int = 50000):
        # Define the data to iterate over
        self._tabix = pysam.TabixFile(source_fn)
        self._current_reader = None  # type: ty.Iterable[ty.Tuple[str, ty.Any, ty.Dict[ty.Any, ty.Any]]]

        # Store the "current" row (sometimes our dbsnp reader will get AHEAD of the sumstats reader,
        #   or be asked for the same position twice in a row, eg multiallelic variants)
        self._current_row = None  # type: ty.Tuple[str, ty.Any, ty.Dict[ty.Any, ty.Any]]

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
        human_chrom_to_tabix = make_chrom_to_contigs(self._tabix)
        dnsnp_chrom = human_chrom_to_tabix[chrom]
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
                self._current_row = next(self._current_reader)  # type: ignore
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
    # Local development (uses part of sumstats, with a matching section of dbsnp)
    # gwas_reader = sniffers.guess_gwas_standard(
    #     '/Users/abought/code/personal/zorp/loaders/benchmarks/data/summary_stats_1k_lines.tab'
    # )
    # rsid_finder = lookups.SnpToRsid(
    #     '/Users/abought/code/personal/zorp/loaders/benchmarks/data/dbsnp_b153_fragment_for_sumstats_benchmarks.lmdb'
    # )
    # rsid_finder = LookupRsidsTabix(
    #     '/Users/abought/code/personal/zorp/loaders/benchmarks/data/dbsnp_b153_fragment_for_sumstats_benchmarks.gz'
    # )

    ###
    # Options useful on the cluster (full dataset)
    gwas_reader = sniffers.guess_gwas_standard('/home/abought/code/zorp/loaders/benchmarks/data/summary_stats.gz')
    # rsid_finder = LookupRsidsTabix('/home/abought/dbsnp/b153/GCF_000001405.25.gz')
    rsid_finder = lookups.SnpToRsid('/home/abought/dbsnp/b153/dbSNP_grch37_b153.lmdb')

    # Perform lookups and track results
    gwas_reader.add_lookup(
        'rsid', lambda variant: rsid_finder(variant.chrom, variant.pos, variant.ref, variant.alt)  # type: ignore
    )

    gwas_reader.write('deleteme.gz')
