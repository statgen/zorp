"""Tests for file format detection"""
import os
import pytest

from zorp import exceptions, readers, sniffers


def _fixture_to_strings(lines: list, delimiter: str = '\t') -> list:
    """Helper so that our unit tests are a little more readable"""
    return [delimiter.join(line)
            for line in lines]


class TestFormatDetection:
    def test_opens_gzip(self):
        fn = os.path.join(os.path.dirname(__file__), "data/sample.gz")
        reader = sniffers.get_reader(fn)
        assert reader is readers.TabixReader

    def test_opens_txt(self):
        fn = os.path.join(os.path.dirname(__file__), "data/pheweb-samples/has-fields-.txt")
        reader = sniffers.get_reader(fn)
        assert reader is readers.TextFileReader

    def test_opens_iterable_if_not_given_a_string(self):
        reader = sniffers.get_reader([['Header', 'Fields'], [1,2]])
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
        assert n == 0, 'Skipped two header rows'
        assert content is None, 'Found correct header row'

    def test_stops_header_search_after_limit(self):
        reader = readers.IterableReader(['walrus', 'carpenter'], parser=None)
        with pytest.raises(exceptions.SnifferException):
            sniffers.get_headers(reader, delimiter='\t', max_check=1)

    def test_no_headers_in_short_file(self):
        reader = readers.IterableReader(['walrus', 'carpenter'], parser=None)
        with pytest.raises(exceptions.SnifferException):
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
    return ('pvalue', 'p.value', 'pval', 'p_score')


class TestFindColumn:
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
        assert actual == {'pval_col': 1, 'is_log_pval': True}

    def test_checks_that_pvalues_are_in_a_realistic_range_0_to_1(self):
        headers = ['pval']
        data = [[100]]
        actual = sniffers.get_pval_column(headers, data)
        assert actual is None


class TestFileFormatDetection:
    def test_warns_if_file_unreadable(self):
        data = _fixture_to_strings([
            ['rsid', 'pval'],
            ['rs1234', '0.5'],
        ])
        with pytest.raises(exceptions.SnifferException):
            sniffers.guess_gwas(data)

    @pytest.mark.skip
    def test_parses_bolt_lmm(self):
        pass

    def test_parses_epacts(self):
        data = _fixture_to_strings([
            ['#CHROM', 'BEGIN', 'END', 'MARKER_ID', 'NS', 'AC', 'CALLRATE', 'MAF', 'PVALUE', 'SCORE', 'N.CASE', 'N.CTRL', 'AF.CASE', 'AF.CTRL'],
            ['20', '1610894', '1610894', '20:1610894_G/A_Synonymous:SIRPG', '266', '138.64', '1', '0.26061', '6.9939e-05', '3.9765', '145', '121', '0.65177', '0.36476']
        ])

        actual = sniffers.guess_gwas(data)
        assert actual._parser._marker_col == 3, 'Found index of marker col'
        assert actual._parser._pval_col == 8, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'

    def test_parses_metal(self):
        data = _fixture_to_strings([
            ['#CHROM', 'POS', 'REF', 'ALT', 'N', 'POOLED_ALT_AF', 'DIRECTION_BY_STUDY', 'EFFECT_SIZE', 'EFFECT_SIZE_SD', 'H2', 'PVALUE'],
            ['1', '10177', 'A', 'AC', '491984', '0.00511094', '?-????????????????-????+???????????????????????????????????????????????????????????????????-????????????????????????????????????????????????????????????????????????????????', '-0.0257947', '0.028959', '1.61266e-06', '0.373073']
        ])
        actual = sniffers.guess_gwas(data)
        assert actual._parser._chr_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'
        assert actual._parser._pval_col == 10, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'

    def test_parses_plink(self):
        # Format: https: // www.cog - genomics.org / plink2 / formats
        # Sample: https: // github.com / babelomics / babelomics / wiki / plink.assoc
        # h/t Josh Weinstock
        data = _fixture_to_strings([
            ['CHR', 'SNP', 'BP', 'A1', 'F_A', 'F_U', 'A2', 'CHISQ', 'P'],
            ['1', 'rs3094315', '742429', 'C', '0.1509', '0.1394', 'T', '0.0759', '0.782', '1.097']
        ])
        actual = sniffers.guess_gwas(data)
        assert actual._parser._chr_col == 0, 'Found index of col'
        assert actual._parser._pos_col == 2, 'Found index of pos col'
        assert actual._parser._ref_col == 3, 'Found index of ref col'
        assert actual._parser._alt_col == 6, 'Found index of alt col'
        assert actual._parser._pval_col == 8, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'

    def test_parses_raremetal(self):
        data = _fixture_to_strings([
            ['#CHROM', 'POS', 'REF', 'ALT', 'N', 'POOLED_ALT_AF', 'DIRECTION_BY_STUDY', 'EFFECT_SIZE', 'EFFECT_SIZE_SD', 'H2', 'PVALUE'],
            ['1', '10177', 'A', 'AC', '491984', '0.00511094', '?-????????????????-????+???????????????????????????????????????????????????????????????????-????????????????????????????????????????????????????????????????????????????????', '-0.0257947', '0.028959', '1.61266e-06', '0.373073']
        ])
        actual = sniffers.guess_gwas(data)
        assert actual._parser._chr_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'
        assert actual._parser._pval_col == 10, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'

    def test_parses_raremetalworker(self):
        data = _fixture_to_strings([
            ['#CHROM', 'POS', 'REF', 'ALT', 'N_INFORMATIVE', 'FOUNDER_AF', 'ALL_AF', 'INFORMATIVE_ALT_AC', 'CALL_RATE', 'HWE_PVALUE', 'N_REF', 'N_HET', 'N_ALT', 'U_STAT', 'SQRT_V_STAT', 'ALT_EFFSIZE', 'PVALUE'],
            ['9', '400066155', 'T', 'C', '432', '0', '0', '0', '1', '1', '432', '0', '0', 'NA', 'NA', 'NA', 'NA']
        ])
        actual = sniffers.guess_gwas(data)
        assert actual._parser._chr_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'
        assert actual._parser._pval_col == 16, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'

    def test_parses_rvtests(self):
        data = _fixture_to_strings([
            ['CHROM', 'POS', 'REF', 'ALT', 'N_INFORMATIVE', 'AF', 'INFORMATIVE_ALT_AC', 'CALL_RATE', 'HWE_PVALUE', 'N_REF', 'N_HET', 'N_ALT', 'U_STAT', 'SQRT_V_STAT', 'ALT_EFFSIZE', 'PVALUE'],
            ['1', '761893', 'G', 'T', '19292', '2.59624e-05:0.000655308:0', '1:1:0', '0.998289:0.996068:0.998381', '1:1:1', '19258:759:18499', '1:1:0', '0:0:0', '1.33113', '0.268484', '18.4664', '7.12493e-07']
        ])
        actual = sniffers.guess_gwas(data)
        assert actual._parser._chr_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 1, 'Found index of pos col'
        assert actual._parser._ref_col == 2, 'Found index of ref col'
        assert actual._parser._alt_col == 3, 'Found index of alt col'
        assert actual._parser._pval_col == 15, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'

    def test_parses_saige(self):
        data = _fixture_to_strings([
            ['CHR', 'POS', 'SNPID', 'Allele1', 'Allele2', 'AC_Allele2', 'AF_Allele2', 'N', 'BETA', 'SE', 'Tstat', 'p.value', 'p.value.NA', 'Is.SPA.converge', 'varT', 'varTstar'],
            ['chr1', '76792', 'chr1:76792:A:C', 'A', 'C', '57', '0.00168639048933983', '16900', '0.573681678183941', '0.663806747906141', '1.30193005902619', '0.387461577915637', '0.387461577915637', '1', '2.2694293866027', '2.41152256615949']
        ])
        actual = sniffers.guess_gwas(data)
        assert actual._parser._marker_col == 2, 'Found index of marker col'
        assert actual._parser._pval_col == 11, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'

    def test_parses_handles_a_mystery_format(self):
        # TODO: Identify the program used and make test more explicit
        # FIXME: This test underscores difficulty of reliable ref/alt detection- a1 comes
        #   before a0, but it might be more valid to switch the order of these columns
        data = _fixture_to_strings([
            ['chr', 'rs', 'ps', 'n_mis', 'n_obs', 'allele1', 'allele0', 'af', 'beta', 'se', 'p_score'],
            ['1', 'rs75333668', '762320', '0', '3610', 'T', 'C', '0.013', '-5.667138e-02', '1.027936e-01', '5.814536e-01']
        ])
        actual = sniffers.guess_gwas(data)
        assert actual._parser._chr_col == 0, 'Found index of chr col'
        assert actual._parser._pos_col == 2, 'Found index of pos col'
        assert actual._parser._ref_col == 5, 'Found index of ref col'
        assert actual._parser._alt_col == 6, 'Found index of alt col'
        assert actual._parser._pval_col == 10, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'

    # Template
    def test_parses_output_of_alisam_pipeline(self):
        data = _fixture_to_strings([
            ['MarkerName', 'chr', 'pos', 'ref', 'alt', 'minor.allele', 'maf', 'mac', 'n', 'pvalue', 'SNPID', 'BETA', 'SE', 'ALTFreq', 'SNPMarker'],
            ['chr1-281876-AC-A', 'chr1', '281876', 'AC', 'A', 'alt', '0.231428578495979', '1053', '2275', '0.447865946615285', 'rs72502741', '-0.0872936159370696', '0.115014743551501', '0.231428578495979', 'chr1:281876_AC/A']
        ])
        actual = sniffers.guess_gwas(data)
        assert actual._parser._chr_col == 1, 'Found index of chr col'
        assert actual._parser._pos_col == 2, 'Found index of pos col'
        assert actual._parser._ref_col == 3, 'Found index of ref col'
        assert actual._parser._alt_col == 4, 'Found index of alt col'
        assert actual._parser._pval_col == 9, 'Found index of pval col'
        assert actual._parser._is_log_pval is False, 'Determined whether is log'
