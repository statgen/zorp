
from collections import abc
import os

import pytest

import parsers
import readers


@pytest.fixture
def simple_tabix_parser():
    return readers.TabixReader(filename=os.path.join(os.path.dirname(__file__), "data/sample.gz"),
                               parser=parsers.NamedFieldsParser(['chr', 'pos', 'ref', 'alt', 'pval']))


@pytest.fixture
def simple_file_parser():
    return readers.FileReader(filename=os.path.join(os.path.dirname(__file__), "data/pheweb-samples/has-fields-.txt"),
                              parser=parsers.NamedFieldsParser(['chr', 'pos', 'ref', 'alt', 'pval']))


def test_tabix_mode_retrieves_data(simple_tabix_parser):
    iterator = simple_tabix_parser.fetch(1, 800_000, 900_000)
    assert isinstance(iterator, abc.Iterable), "returns an iterator"
    all_records = list(iterator)
    assert len(all_records) == 1, "one record in region"

    # Then try a different region. Should reuse same iterator.
    iter2 = simple_tabix_parser.fetch("X", 2_600_000, 3_000_000)
    all_records = list(iter2)
    assert len(all_records) == 3, "one record in region"


def test_filemode_iterator(simple_file_parser):
    assert isinstance(simple_file_parser, abc.Iterable)

    iterator = iter(simple_file_parser)
    first_row = next(iterator)
    assert isinstance(first_row, tuple), "rows are parsed as tuples"


@pytest.mark.skip
def test_filter_criteria_with_tabix(simple_tabix_parser):
    # assert False
    pass


def test_filter_criteria_with_method(simple_file_parser):
    # Can define a filter using any function
    simple_file_parser.add_filter("chr", lambda val, row: val == "1")
    assert len(simple_file_parser._filters) == 1

    # File will act on it
    assert len(list(simple_file_parser)) == 7, "output was restricted to the expected rows"


def test_filter_criteria_with_value(simple_file_parser):
    # Can define a filter that exactly matches a value
    simple_file_parser.add_filter("chr", "1")
    assert len(simple_file_parser._filters) == 1

    # File will act on it
    assert len(list(simple_file_parser)) == 7, "output was restricted to the expected rows"


def test_must_specify_filename_or_source():
    with pytest.raises(Exception):
        readers.FileReader(filename="anything", source=[])


def test_can_specify_iterable_as_source():
    # Any row-based iterable can serve as the source and will be parsed as expected
    parser = readers.FileReader(source=["X\t1\tA\tG"])
    result = next(iter(parser))
    assert result[0] == "X"


def test_sorting_stuff():
    parser = readers.FileReader(filename=os.path.join(os.path.dirname(__file__), "data/unsorted.txt"))
    prev_chr = 0
    prev_pos = 0
    for row in parser.sort(1, 2):
        assert int(row[0]) > prev_chr, "Chromosomes always increase"
        assert int(row[1]) > prev_pos, "Positions always increase"


# TODO: Expected api: [ row['chr'] for row in source(filename).sort().add_filter("chr", lambda val: value==1) ]