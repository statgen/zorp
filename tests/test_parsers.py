"""Test line parsers"""

import pytest

from zorp import (
    exceptions,
    parsers,
)


class TestTupleParser:
    def test_can_configure_delimiter(self):
        parser = parsers.TupleLineParser(delimiter=',')
        output = parser('a,b,c')
        assert output == ('a', 'b', 'c')

    def test_can_return_list_instead_of_tuple(self):
        parser = parsers.TupleLineParser(container=list)
        output = parser('a\tb\tc')
        assert isinstance(output, list)
        assert output == ['a', 'b', 'c']


class TestStandardGwasParser:
    def test_parses_locuszoom_standard_format(self):
        line = '1\t100\tA\tC\t10'
        output = parsers.standard_gwas_parser(line)
        assert isinstance(output, parsers._basic_standard_container)
        assert output.chrom == '1'
        assert output.pos == 100
        assert output.ref == 'A'
        assert output.alt == 'C'
        assert output.pvalue == pytest.approx(1e-10)

    def test_enforces_pos_as_int(self):
        line = '1\tNOPE\tA\tC\t0.05'
        with pytest.raises(exceptions.LineParseException, match="invalid literal"):
            parsers.standard_gwas_parser(line)

    def test_enforces_readable_pvalue(self):
        line = '1\t100\tA\tC\tNOPE'
        with pytest.raises(exceptions.LineParseException, match="could not convert string to float"):
            parsers.standard_gwas_parser(line)

    def test_handles_missing_pvalues(self):
        line = '1\t100\tA\tC\tNA'
        p = parsers.standard_gwas_parser(line)
        assert p.pvalue is None, "Fills placeholder values with python None"

    def test_can_convert_to_log(self):
        line = '1\t100\tA\tC\t1'
        special_parser = parsers.GenericGwasLineParser(chr_col=1, pos_col=2, ref_col=3, alt_col=4,
                                                       pval_col=5, is_log_pval=True,
                                                       delimiter='\t')
        p = special_parser(line)
        assert p.log_pvalue == pytest.approx(1), 'Converts -log to pvalue'
        assert p.pvalue == pytest.approx(0.1), 'Converts -log to pvalue'

    def test_parses_marker_to_clean_format(self):
        line = 'chr2:100:A:C_anno\t.05'
        special_parser = parsers.GenericGwasLineParser(marker_col=1, pval_col=2, delimiter='\t')
        p = special_parser(line)
        assert p.chrom == '2', 'Finds chromosome'
        assert p.pos == 100, 'Finds position'
        assert p.ref == 'A', 'Finds ref'
        assert p.alt == 'C', 'Finds alt'
        assert p.marker == '2:100_A/C', 'Turns a messy marker into a cleaned standardized format'
