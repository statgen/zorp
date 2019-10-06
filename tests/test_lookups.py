"""Tests for lookup functionality"""

import os

import pytest

from zorp import lookups


@pytest.fixture
def rsid_testdata():
    data_fn = os.path.join(os.path.dirname(__file__), "data/rsid-build37-segment.lmdb")
    return lookups.FindRsid(data_fn)


class TestFindRsid:
    def test_finds_snp_in_lookup(self, rsid_testdata):
        res = rsid_testdata('1', 278_829, 'T', 'G')
        assert res == 'rs1415452960'

    # def test_finds_several_snps_in_sequence(self, rsid_testdata):
    #     # Ensures that no weird state prevents the closure from working > 1x
    #     # TODO: this will allow us to make optimizations later
    #     scenarios = [
    #
    #     ]
    #     for args, rs in scenarios:
    #         assert rsid_testdata(*args) == rs

    def test_distinguishes_between_two_rsids_same_position(self, rsid_testdata):
        res = rsid_testdata('1', 10_051, 'A', 'AC')
        assert res == 'rs1326880612'

        res = rsid_testdata('1', 10_051, 'A', 'G')
        assert res == 'rs1052373574'

    def test_generation_script_doesnt_drop_last_row(self, rsid_testdata):
        # Spot check that the "consolidate rows" mechanism doesn't lose last row
        res = rsid_testdata('1', 600349, 'A', 'G')
        assert res == 'rs1006346155'
