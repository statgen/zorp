"""Tests for lookup functionality"""

import os

import pytest

from zorp import lookups


@pytest.fixture
def rsid_testdata():
    data_fn = os.path.join(os.path.dirname(__file__), "data/rsid_lookup/dbSNP_grch38_b153.lmdb")
    return lookups.FindRsid(data_fn)


class TestFindRsid:
    def test_finds_snp_in_lookup(self, rsid_testdata):
        res = rsid_testdata('1', 10063, 'A', 'C')
        assert res == 1010989343

    def test_identifies_two_alts_with_same_rsid(self, rsid_testdata):
        res = rsid_testdata('1', 10164, 'A', 'G')
        assert res == 1413947121

        res = rsid_testdata('1', 10164, 'A', 'T')
        assert res == 1413947121

    def test_distinguishes_two_rsids_for_same_position(self, rsid_testdata):
        res = rsid_testdata('1', 10051, 'A', 'G')
        assert res == 1052373574

        res = rsid_testdata('1', 10051, 'A', 'AC')
        assert res == 1326880612
