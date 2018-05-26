
from collections import abc
import collections
import os

import pytest

import parsers


class SimpleFileParser(parsers.BaseGWASParser):
    __one_row__ = collections.namedtuple("gwas_data", ["chr", "pos", "ref", "alt", "pvalue"])

    def _parse_row(self, row: str):
        row = row.split("\t")
        return self.__one_row__(row[0], None, None, None, None)


@pytest.fixture
def simple_tabix():
    return SimpleFileParser(filename=os.path.join(os.path.dirname(__file__), "data/sample.gz"))


@pytest.fixture
def simple_file():
    return SimpleFileParser(filename=os.path.join(os.path.dirname(__file__), "data/pheweb-samples/has-fields-.txt"))


def test_tabix_mode_retrieves_data(simple_tabix):
    iterator = simple_tabix.fetch(1, 800_000, 900_000)
    assert isinstance(iterator, abc.Iterable), "returns an iterator"
    all_records = list(iterator)
    assert len(all_records) == 1, "one record in region"

    # Then try a different region. Should reuse same iterator.
    iter2 = simple_tabix.fetch("X", 2_600_000, 3_000_000)
    all_records = list(iter2)
    assert len(all_records) == 3, "one record in region"


def test_filemode_iterator(simple_file):
    assert isinstance(simple_file, abc.Iterable)

    iterator = iter(simple_file)
    first_row = next(iterator)
    assert isinstance(first_row, tuple), "rows are parsed as tuples"


def test_filter_criteria_with_tabix(simple_tabix):
    # assert False
    pass


def test_filter_criteria_with_file(simple_file):
    # Can define a filter
    simple_file.add_filter("chr", lambda val,row: val == "1")
    assert len(simple_file._filters) == 1

    # File will act on it
    assert len(list(simple_file)) == 7, "output was restricted to the expected rows"


def test_must_specify_filename_or_source():
    with pytest.raises(Exception):
        SimpleFileParser(filename="anything", source=[])


def test_can_specify_iterable_as_source():
    # Any row-based iterable can serve as the source and will be parsed as expected
    parser = SimpleFileParser(source=["X\t1\tA\tG"])
    result = next(iter(parser))
    assert result[0] == "X"


# TODO: Expected api: [ row['chr'] for row in source(filename).sort().add_filter("chr", lambda val: value==1) ]