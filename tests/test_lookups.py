"""Tests for lookup functionality"""

import os

import pytest

from zorp import lookups


@pytest.fixture
def rsid_testdata():
    data_fn = os.path.join(os.path.dirname(__file__), "data/snp_to_rsid/dbSNP_grch38_b153.lmdb")
    return lookups.SnpToRsid(data_fn)


class TestSnpToRsid:
    def test_finds_snp_in_lookup(self, rsid_testdata):
        res = rsid_testdata('1', 10063, 'A', 'C')
        assert res == 'rs1010989343'

    def test_finds_no_snp_if_wrong_alt(self, rsid_testdata):
        res = rsid_testdata('1', 10063, 'A', 'A')
        assert res is None

    def test_identifies_two_alts_with_same_rsid(self, rsid_testdata):
        res = rsid_testdata('1', 10164, 'A', 'G')
        assert res == 'rs1413947121'

        res = rsid_testdata('1', 10164, 'A', 'T')
        assert res == 'rs1413947121'

    def test_distinguishes_two_rsids_for_same_position(self, rsid_testdata):
        res = rsid_testdata('1', 10051, 'A', 'G')
        assert res == 'rs1052373574'

        res = rsid_testdata('1', 10051, 'A', 'AC')
        assert res == 'rs1326880612'

    def test_helper_lists_all_known_chroms(self, rsid_testdata):
        known = rsid_testdata.known_chroms()
        assert known == sorted([  # The full set used to build the test data
            '1', '10',
            '11', '12', '13', '14', '15', '16', '17', '18', '19',
            '2', '3', '4', '5', '6', '7', '8', '9', '20',
            '21', '22', 'MT', 'X', 'Y',
        ])
