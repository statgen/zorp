"""Test line parsers"""

import pytest

import exceptions
import parsers


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
        line = '1\t100\tA\tC\t0.05'
        output = parsers.standard_gwas_parser(line)
        assert isinstance(output, parsers._basic_standard_container)
        assert output.chrom == '1'
        assert output.pos == 100
        assert output.ref == 'A'
        assert output.alt == 'C'
        assert output.pvalue == 0.8912509381337456  # standard format provides logpvalues!

    def test_enforces_pos_as_int(self):
        line = '1\tNOPE\tA\tC\t0.05'
        with pytest.raises(exceptions.LineParseException, match="invalid literal"):
            output = parsers.standard_gwas_parser(line)

    def test_enforces_readable_pvalue(self):
        line = '1\t100\tA\tC\tNOPE'
        with pytest.raises(exceptions.LineParseException, match="could not convert string to float"):
            output = parsers.standard_gwas_parser(line)

# TODO: Add tests for the generic gwas parser (and all its various combinations of options)
