"""Tests for file format detection"""

from zorp import sniffers


class TestHeaders:
    def test_scenarios(self):
        assert sniffers.is_header('#Comment') is True, 'Comment lines are headers!'
        assert sniffers.is_header('Header\tLabels') is True, 'Headers tend to be text'
        assert sniffers.is_header('X\t100') is False, 'Data has numbers'
        assert sniffers.is_header('X\t.') is False, 'Missing data is still data'

        assert sniffers.is_header('X,100', delimiter=',') is False, 'Handles data as csv'
        assert sniffers.is_header('//100', comment_char='//') is True, 'Handles different comments'
