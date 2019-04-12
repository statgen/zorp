"""Tests for file format detection"""
import pytest

from zorp import exceptions, sniffers


def _fixture_to_strings(lines: list, delimiter: str ='\t') -> list:
    """Helper so that our unit tests are a little more readable"""
    return [delimiter.join(line)
            for line in lines]


def test_header_detection():
    assert sniffers.is_header('#Comment') is True, 'Comment lines are headers!'
    assert sniffers.is_header('Header\tLabels') is True, 'Headers tend to be text'
    assert sniffers.is_header('X\t100') is False, 'Data has numbers'
    assert sniffers.is_header('X\t.') is False, 'Missing data is still data'

    assert sniffers.is_header('X,100', delimiter=',') is False, 'Handles data as csv'
    assert sniffers.is_header('//100', comment_char='//') is True, 'Handles different comments'


def test_is_numeric():
    assert sniffers.is_numeric('1.23') is True, 'Number is numeric'
    assert sniffers.is_numeric('inf') is True, 'Infinity is numeric'
    assert sniffers.is_numeric('1.23e-4') is True, 'Exponent is numeric'
    assert sniffers.is_numeric('NA') is True, 'Missing value counts as numeric'
    assert sniffers.is_numeric('Antarctica') is False, 'String is not numeric'
    assert sniffers.is_numeric('1.23.4') is False, 'Version string is not numeric'


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
        pass
