"""
Look up the RSID associated with a variant, based on a tabixed copy of dbSNP
"""
import gzip
import typing as ty

from zorp import parsers, readers

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

    def __call__(self, chrom, pos, ref, alt) -> ty.Union[int, None]:
        # Look up the rsid associated with a current SNP.
        self._advance_current_reader(chrom, pos)

        snp_chrom, snp_pos, ref_alt_options = self._current_row
        if snp_chrom != chrom or snp_pos != pos:
            # Ensure that the dbSNP reader is at the correct position!
            return None
        return ref_alt_options.get('{}/{}'.format(ref, alt), None)

    def __del__(self):
        self._vcf.close()


if __name__ == '__main__':
    standard_gwas_parser = parsers.GenericGwasLineParser(
        chrom_col=1, pos_col=2, ref_col=3, alt_col=4,
        pvalue_col=5, is_neg_log_pvalue=True,
        beta_col=6, stderr_beta_col=7,
        allele_freq_col=8,
        is_alt_effect=True,
        delimiter='\t'
    )
    gwas_reader = readers.TabixReader('data/summary_stats.gz', parser=standard_gwas_parser)
    rsid_finder = LookupRsidsTabix('/home/abought/dbsnp/b153/GCF_000001405.25.gz')

    # Perform lookups and track results
    gwas_reader.add_lookup('rsid', rsid_finder)
    gwas_reader.write('deleteme.gz')
