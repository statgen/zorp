"""Tests for file format detection"""
import os
import pytest

from zorp import exceptions, parsers, readers, sniffers


def _fixture_to_strings(lines: list, delimiter: str = '\t') -> list:
    """
    Helper so that our unit tests are a little more readable. Real tabix files give delimited strings,
        not lists of fields
    """
    return [delimiter.join(line)
            for line in lines]


class TestFiletypeDetection:
    def test_opens_gzip(self):
        fn = os.path.join(os.path.dirname(__file__), "data/sample.gz")
        reader = sniffers.get_reader(fn)
        assert reader is readers.TabixReader

    def test_opens_txt(self):
        fn = os.path.join(os.path.dirname(__file__), "data/pheweb-samples/has-fields-.txt")
        reader = sniffers.get_reader(fn)
        assert reader is readers.TextFileReader

    def test_opens_iterable_if_not_given_a_string(self):
        reader = sniffers.get_reader([['Header', 'Fields'], [1, 2]])
        assert reader is readers.IterableReader


class TestHeaderDetection:
    def test_header_detection(self):
        assert sniffers.is_header('#Comment') is True, 'Comment lines are headers!'
        assert sniffers.is_header('Header\tLabels') is True, 'Headers tend to be text'
        assert sniffers.is_header('X\t100') is False, 'Data has numbers'
        assert sniffers.is_header('X\t.') is False, 'Missing data is still data'

        assert sniffers.is_header('X,100', delimiter=',') is False, 'Handles data as csv'
        assert sniffers.is_header('//100', comment_char='//') is True, 'Handles different comments'

    #####
    # Convenience method: Automatic header detection
    def test_can_find_headers(self):
        reader = readers.IterableReader(["#Comment line", '#Header\tLabels', 'X\t100'], parser=None)
        n, content = sniffers.get_headers(reader, delimiter='\t')
        assert n == 2, 'Skipped two header rows'
        assert content == '#Header\tLabels', 'Found correct header row'

    def test_handles_lack_of_headers(self):
        reader = readers.IterableReader(['X\t100', 'X\t101'], parser=None)
        n, content = sniffers.get_headers(reader, delimiter='\t')
        assert n == 0, 'File has no header rows'
        assert content is None, 'No header row, so headers are blank'

    def test_stops_header_search_after_limit(self):
        reader = readers.IterableReader(['walrus', 'carpenter'], parser=None)
        with pytest.raises(exceptions.SnifferException, match='after limit'):
            sniffers.get_headers(reader, delimiter='\t', max_check=1)

    def test_no_headers_in_short_file(self):
        reader = readers.IterableReader(['walrus', 'carpenter'], parser=None)
        with pytest.raises(exceptions.SnifferException, match='entire file'):
            sniffers.get_headers(reader, delimiter='\t')


def test_is_numeric():
    assert sniffers.is_numeric('1.23') is True, 'Number is numeric'
    assert sniffers.is_numeric('inf') is True, 'Infinity is numeric'
    assert sniffers.is_numeric('1.23e-4') is True, 'Exponent is numeric'
    assert sniffers.is_numeric('NA') is True, 'Missing value counts as numeric'
    assert sniffers.is_numeric('Antarctica') is False, 'String is not numeric'
    assert sniffers.is_numeric('1.23.4') is False, 'Version string is not numeric'


@pytest.fixture
def pval_names():
    return 'pvalue', 'p.value', 'pval', 'p_score'


class TestFindColumn:
    def test_levenshtein_handles_empty_string(self):
        assert sniffers.levenshtein('', 'bob') == 3, 'Calculates differences if one string is empty'

    def test_finds_first_exact_match_for_synonym(self, pval_names):
        headers = ['chr', 'pos', 'p.value', 'marker']
        match = sniffers.find_column(pval_names, headers)
        assert match == 2

    def test_chooses_first_exact_match_when_more_than_one_is_present(self, pval_names):
        headers = ['chr', 'pvalue', 'p.value', 'marker']
        match = sniffers.find_column(pval_names, headers)
        assert match == 1

    def test_prefers_exact_matches_over_fuzzy_matches(self, pval_names):
        headers = ['chr1', 'pos1', 'pvalues', 'p.value', '1marker']
        match = sniffers.find_column(pval_names, headers)
        assert match == 3

    def test_finds_the_first_header_that_closely_matches_a_synonym(self, pval_names):
        headers = ['chr', 'pos', 'marker', 'p-value']
        match = sniffers.find_column(pval_names, headers)
        assert match == 3

    def test_returns_none_if_no_good_match_can_be_found(self, pval_names):
        headers = ['chr', 'pos', 'marker', 'pval_score']
        match = sniffers.find_column(pval_names, headers)
        assert match is None

    def test_will_match_based_on_a_configurable_threshold(self, pval_names):
        headers = ['chr', 'marker', 'pval_score']
        match = sniffers.find_column(pval_names, headers, threshold=3)
        assert match == 2

    def test_skips_headers_with_a_null_value(self, pval_names):
        headers = ['chr', None, 'marker', 'pval']
        match = sniffers.find_column(pval_names, headers)
        assert match == 3


class TestGetPvalColumn:
    def test_finds_logp_before_p(self):
        headers = ['logpvalue', 'pval']
        data = [[0.5, 0.5]]
        actual = sniffers.get_pval_column(headers, data)
        assert actual == {'pvalue_col': 1, 'is_neg_log_pvalue': True}

    def test_checks_that_pvalues_are_in_a_realistic_range_0_to_1(self):
        headers = ['pval']
        data = [[100]]
        actual = sniffers.get_pval_column(headers, data)
        assert actual is None


class TestGetEffectSizeColumns:
    def test_invalid_data_doesnt_count_as_effect_size(self):
        headers = ['beta']
        data = [['bork']]
        actual = sniffers.get_effect_size_columns(headers, data)
        assert actual is None


class TestFileFormatDetection:
    def test_warns_if_file_lacks_required_fields(self):
        data = _fixture_to_strings([
            ['rsid', 'pval'],
            ['rs1234', '0.5'],
        ])
        with pytest.raises(exceptions.SnifferException):
            sniffers.guess_gwas_generic(data)

    def test_sniffer_warns_if_cant_find_pval(self):
        data = _fixture_to_strings([
            ['rsid', 'marker'],
            ['rs1234', '0.5'],
        ])
        with pytest.raises(exceptions.SnifferException, match='pvalue'):
            sniffers.guess_gwas_generic(data)

    def test_can_provide_extra_options_for_parser(self):
        data = _fixture_to_strings([
            ['#chrom', 'pos', 'ref', 'alt', 'neg_log_pvalue', 'alt_allele_freq'],
            ['1', '762320', 'C', 'T', '0.36947042857317597', '0.5', '0.1']
        ])
        actual = sniffers.guess_gwas_generic(data, parser_options={'allele_freq_col': 6})
        assert actual._parser._allele_freq_col == 5, 'Sniffer used an option that it could not have auto-detected'

    def test_sniffer_validates_options(self):
        with pytest.raises(exceptions.ConfigurationException, match='exclusive'):
            sniffers.guess_gwas_generic(['1', '2'],
                                        parser=parsers.standard_gwas_parser_basic,
                                        parser_options={'option': 1})

    # Sample file formats/ sample data that Zorp should be able to detect
    def test_can_guess_standard_format(self):
        # Tracks the "standard format" defined as a convenience parser
        data = _fixture_to_strings([
            ['#chrom', 'pos', 'ref', 'alt', 'neg_log_pvalue', 'beta', 'stderr_beta', 'alt_allele_freq'],
            ['1', '762320', 'C', 'T', '0.36947042857317597', '0.5', '0.1', '0.5']
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._chrom_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'

        assert actual._parser._pvalue_col == 4, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is True, 'Determined whether is log'

        assert actual._parser._beta_col == 5, 'beta field detected'
        assert actual._parser._stderr_col == 6, 'stderr_beta field detected'

    def test_can_guess_bolt_lmm(self):
        data = _fixture_to_strings([
            ['SNP', 'CHR', 'BP', 'A1', 'A0', 'MAF', 'HWEP', 'INFO', 'BETA', 'SE', 'P'],
            ['10:48698435_A_G', '10', '48698435', 'A', 'G', '0.01353', '0.02719', '0.960443', '0.0959329', '0.0941266', '3.3E-01']  # noqa: E501
        ])

        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._marker_col == 0, 'Found index of marker col'
        assert actual._parser._pvalue_col == 10, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 8, 'beta field detected'
        assert actual._parser._stderr_col == 9, 'stderr_beta field detected'

    def test_can_guess_epacts(self):
        data = _fixture_to_strings([
            ['#CHROM', 'BEGIN', 'END', 'MARKER_ID', 'NS', 'AC', 'CALLRATE', 'MAF', 'PVALUE', 'SCORE', 'N.CASE', 'N.CTRL', 'AF.CASE', 'AF.CTRL'],  # noqa: E501
            ['20', '1610894', '1610894', '20:1610894_G/A_Synonymous:SIRPG', '266', '138.64', '1', '0.26061', '6.9939e-05', '3.9765', '145', '121', '0.65177', '0.36476']  # noqa: E501
        ])

        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._marker_col == 3, 'Found index of marker col'
        assert actual._parser._pvalue_col == 8, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col is None, 'No beta field detected'
        assert actual._parser._stderr_col is None, 'No stderr_beta field detected'

    def can_guess_emmax_epacts(self):
        """Fileformat sample provided by multiple tools"""
        data = _fixture_to_strings([
            ['#CHROM', 'BEG', 'END', 'MARKER_ID', 'NS', 'AC', 'CALLRATE', 'GENOCNT', 'MAF', 'STAT', 'PVALUE', 'BETA', 'SEBETA', 'R2'],  # noqa: E501
            ['1', '762320', '762320', '1:762320_C/T_rs75333668', '3805', '100.00', '1.00000', '3707/96/2', '0.01314', '0.7942', '0.4271', '0.08034', '0.1012', '0.0001658']  # noqa: E501
        ])

        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._marker_col == 3, 'Found index of marker col'
        assert actual._parser._pvalue_col == 10, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 11, 'beta field detected'
        assert actual._parser._stderr_col == 12, 'stderr_beta field detected'

    def test_can_guess_metal(self):
        data = _fixture_to_strings([
            ['#CHROM', 'POS', 'REF', 'ALT', 'N', 'POOLED_ALT_AF', 'DIRECTION_BY_STUDY', 'EFFECT_SIZE', 'EFFECT_SIZE_SD', 'H2', 'PVALUE'],  # noqa: E501
            ['1', '10177', 'A', 'AC', '491984', '0.00511094', '?-????????????????-????+???????????????????????????????????????????????????????????????????-????????????????????????????????????????????????????????????????????????????????', '-0.0257947', '0.028959', '1.61266e-06', '0.373073']  # noqa: E501
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._chrom_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'
        assert actual._parser._pvalue_col == 10, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 7, 'beta field detected'
        assert actual._parser._stderr_col == 8, 'stderr_beta field detected'

    def test_can_guess_plink(self):
        # Format: https://www.cog-genomics.org/plink2/formats
        # Sample: https://github.com/babelomics/babelomics/wiki/plink.assoc
        # h/t Josh Weinstock
        data = _fixture_to_strings([
            ['CHR', 'SNP', 'BP', 'A1', 'F_A', 'F_U', 'A2', 'CHISQ', 'P'],
            ['1', 'rs3094315', '742429', 'C', '0.1509', '0.1394', 'T', '0.0759', '0.782', '1.097']
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._chrom_col == 0, 'Found index of col'
        assert actual._parser._pos_col == 2, 'Found index of pos col'
        assert actual._parser._ref_col == 3, 'Found index of ref col'
        assert actual._parser._alt_col == 6, 'Found index of alt col'
        assert actual._parser._pvalue_col == 8, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col is None, 'No beta field detected'
        assert actual._parser._stderr_col is None, 'No stderr_beta field detected'

    def test_can_guess_raremetal(self):
        data = _fixture_to_strings([
            ['#CHROM', 'POS', 'REF', 'ALT', 'N', 'POOLED_ALT_AF', 'DIRECTION_BY_STUDY', 'EFFECT_SIZE', 'EFFECT_SIZE_SD', 'H2', 'PVALUE'],  # noqa: E501
            ['1', '10177', 'A', 'AC', '491984', '0.00511094', '?-????????????????-????+???????????????????????????????????????????????????????????????????-????????????????????????????????????????????????????????????????????????????????', '-0.0257947', '0.028959', '1.61266e-06', '0.373073']  # noqa: E501
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._chrom_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'
        assert actual._parser._pvalue_col == 10, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 7, 'Beta field detected'
        assert actual._parser._stderr_col == 8, 'stderr_beta field detected'

    def test_can_guess_raremetalworker(self):
        data = _fixture_to_strings([
            ['#CHROM', 'POS', 'REF', 'ALT', 'N_INFORMATIVE', 'FOUNDER_AF', 'ALL_AF', 'INFORMATIVE_ALT_AC', 'CALL_RATE', 'HWE_PVALUE', 'N_REF', 'N_HET', 'N_ALT', 'U_STAT', 'SQRT_V_STAT', 'ALT_EFFSIZE', 'PVALUE'],  # noqa: E501
            ['9', '400066155', 'T', 'C', '432', '0', '0', '0', '1', '1', '432', '0', '0', 'NA', 'NA', 'NA', 'NA']
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._chrom_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'
        assert actual._parser._pvalue_col == 16, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 15, 'beta field detected'
        assert actual._parser._stderr_col is None, 'No stderr_beta field detected'

    def test_can_guess_rvtests(self):
        data = _fixture_to_strings([
            ['CHROM', 'POS', 'REF', 'ALT', 'N_INFORMATIVE', 'AF', 'INFORMATIVE_ALT_AC', 'CALL_RATE', 'HWE_PVALUE', 'N_REF', 'N_HET', 'N_ALT', 'U_STAT', 'SQRT_V_STAT', 'ALT_EFFSIZE', 'PVALUE'],  # noqa: E501
            ['1', '761893', 'G', 'T', '19292', '2.59624e-05:0.000655308:0', '1:1:0', '0.998289:0.996068:0.998381', '1:1:1', '19258:759:18499', '1:1:0', '0:0:0', '1.33113', '0.268484', '18.4664', '7.12493e-07']  # noqa: E501
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._chrom_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'
        assert actual._parser._pvalue_col == 15, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 14, 'beta field detected'
        assert actual._parser._stderr_col is None, 'No stderr_beta field detected'

    def test_can_guess_saige(self):
        data = _fixture_to_strings([
            ['CHR', 'POS', 'SNPID', 'Allele1', 'Allele2', 'AC_Allele2', 'AF_Allele2', 'N', 'BETA', 'SE', 'Tstat', 'p.value', 'p.value.NA', 'Is.SPA.converge', 'varT', 'varTstar'],  # noqa: E501
            ['chr1', '76792', 'chr1:76792:A:C', 'A', 'C', '57', '0.00168639048933983', '16900', '0.573681678183941', '0.663806747906141', '1.30193005902619', '0.387461577915637', '0.387461577915637', '1', '2.2694293866027', '2.41152256615949']  # noqa: E501
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._marker_col == 2, 'Found index of marker col'
        assert actual._parser._pvalue_col == 11, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 8, 'beta field detected'
        assert actual._parser._stderr_col == 9, 'stderr_beta field detected'

    def test_can_guess_a_mystery_format(self):
        # TODO: Identify the program used and make test more explicit
        # FIXME: This test underscores difficulty of reliable ref/alt detection- a1 comes
        #   before a0, but it might be more valid to switch the order of these columns. Leave meaning up to the user.
        data = _fixture_to_strings([
            ['chr', 'rs', 'ps', 'n_mis', 'n_obs', 'allele1', 'allele0', 'af', 'beta', 'se', 'p_score'],
            ['1', 'rs75333668', '762320', '0', '3610', 'T', 'C', '0.013', '-5.667138e-02', '1.027936e-01', '5.814536e-01']  # noqa: E501
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._chrom_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 2, 'Found index of pos col'
        assert actual._parser._ref_col == 5, 'Found index of ref col'
        assert actual._parser._alt_col == 6, 'Found index of alt col'
        assert actual._parser._pvalue_col == 10, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 8, 'beta field detected'
        assert actual._parser._stderr_col == 9, 'stderr_beta field detected'

    def test_can_guess_output_of_alisam_pipeline(self):
        data = _fixture_to_strings([
            ['MarkerName', 'chr', 'pos', 'ref', 'alt', 'minor.allele', 'maf', 'mac', 'n', 'pvalue', 'SNPID', 'BETA', 'SE', 'ALTFreq', 'SNPMarker'],  # noqa: E501
            ['chr1-281876-AC-A', 'chr1', '281876', 'AC', 'A', 'alt', '0.231428578495979', '1053', '2275', '0.447865946615285', 'rs72502741', '-0.0872936159370696', '0.115014743551501', '0.231428578495979', 'chr1:281876_AC/A']  # noqa: E501
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._chrom_col == 1, 'Found index of chr col'
        assert actual._parser._pos_col == 2, 'Found index of pos col'
        assert actual._parser._ref_col == 3, 'Found index of ref col'
        assert actual._parser._alt_col == 4, 'Found index of alt col'
        assert actual._parser._pvalue_col == 9, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 11, 'beta field detected'
        assert actual._parser._stderr_col == 12, 'stderr_beta field detected'

    def test_can_guess_whatever_diagram_was_using(self):
        # FIXME: If this format turns out to be common, we should improve it to fetch all four values, instead of just
        #   the two that the marker will provide
        data = _fixture_to_strings([
            ['Chr:Position', 'Allele1', 'Allele2', 'Effect', 'StdErr', 'P-value', 'TotalSampleSize'],
            ['5:29439275', 'T', 'C', '-0.0003', '0.015', '0.99', '111309'],
        ])
        actual = sniffers.guess_gwas_generic(data)
        assert actual._parser._marker_col == 0, 'Found index of marker col'
        assert actual._parser._pvalue_col == 5, 'Found index of pval col'
        assert actual._parser._is_neg_log_pvalue is False, 'Determined whether is log'

        assert actual._parser._beta_col == 3, 'beta field detected'
        assert actual._parser._stderr_col == 4, 'stderr_beta field detected'
