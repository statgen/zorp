import pytest

from zorp import parsers


##
# Sample preconfigured parsers that work with the data in these test suites
# These were defined at a time before optional columns were introduced, so in practice many of our tests have the
#   same layout and order of fields
@pytest.fixture(scope='module')
def standard_gwas_parser():
    return parsers.GenericGwasLineParser(chrom_col=1, pos_col=2, ref_col=3, alt_col=4,
                                         pvalue_col=5, is_neg_log_pvalue=True,
                                         beta_col=6, stderr_beta_col=7,
                                         allele_freq_col=8,
                                         is_alt_effect=True,
                                         delimiter='\t')


@pytest.fixture(scope='module')
def standard_gwas_parser_basic():
    return parsers.GenericGwasLineParser(chrom_col=1, pos_col=2, ref_col=3, alt_col=4,
                                         pvalue_col=5, is_neg_log_pvalue=True,
                                         delimiter='\t')
