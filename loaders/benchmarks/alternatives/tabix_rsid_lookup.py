#!/usr/bin/env python
"""
Look up the RSID associated with a variant, based on a tabixed copy of dbSNP
"""
import gzip
import typing as ty

from zorp import parsers, sniffers

from loaders.make_rsid_lookup import line_parser, make_group_iterator, make_file_iterator


class LookupRsidsTabix:
    """
    Find RSIDs that match a specified file. This is a tabix-based design and inherently assumes that each lookup
        performed will be sequential (after the previous one)
    """
    def __init__(self, source_fn: str):
        # Define the data to iterate over
        self._vcf = gzip.open(source_fn, "rt")
        all_lines = make_file_iterator(self._vcf)
        self._position_entries = iter(make_group_iterator(all_lines))

        # Store the "current" row (sometimes our dbsnp reader will get AHEAD of the sumstats reader,
        #   or be asked for the same position twice in a row, eg multiallelic variants)
        self._current_row = next(self._position_entries)

    def _advance_current_reader(self, target_chrom, target_pos):
        chrom, pos, _ = self._current_row

        if chrom == target_chrom and pos >= target_pos:
            # Don't advance the dbSNP reader any further- we are already at or past where we need to be
            return self._current_row

        # If the dbSNP reader is lagging the gwas reader, then we need to advance
        while chrom != target_chrom or pos < target_pos:
            try:
                self._current_row = next(self._position_entries)
            except StopIteration:
                # If we hit the end of the dbSNP file, exit the loop and never let the iterator advance again
                break

            chrom, pos, _ = self._current_row

    def __call__(self, variant: parsers.BasicVariant) -> ty.Union[int, None]:
        # Look up the rsid associated with a current SNP.
        self._advance_current_reader(variant.chrom, variant.pos)

        snp_chrom, snp_pos, ref_alt_options = self._current_row
        if snp_chrom != variant.chrom or snp_pos != variant.pos:
            # Ensure that the dbSNP reader is at the correct position!
            return None
        return ref_alt_options.get('{}/{}'.format(variant.ref, variant.alt), None)

    def __del__(self):
        if self._vcf:
            self._vcf.close()


if __name__ == '__main__':
    gwas_reader = sniffers.guess_gwas_standard('/home/abought/code/zorp/loaders/benchmarks/data/summary_stats.gz')
    rsid_finder = LookupRsidsTabix('/home/abought/dbsnp/b153/GCF_000001405.25.gz')

    # Perform lookups and track results
    gwas_reader.add_lookup('rsid', rsid_finder)
    gwas_reader.write('deleteme.gz')
