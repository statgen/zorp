"""
Test reader functionality
"""
from collections import abc
import os

import pytest

import parsers
import readers


@pytest.fixture
def simple_tabix_reader():
    return readers.TabixReader(filename=os.path.join(os.path.dirname(__file__), "data/sample.gz"),
                               parser=parsers.standard_gwas_parser,
                               skip_rows=1)


@pytest.fixture
def simple_file_reader():
    return readers.TextFileReader(filename=os.path.join(os.path.dirname(__file__), "data/pheweb-samples/has-fields-.txt"),
                                  parser=parsers.standard_gwas_parser,
                                  skip_rows=1)


class TestTabixReader:
    def test_tabix_mode_retrieves_data(self, simple_tabix_reader):
        iterator = simple_tabix_reader.fetch(1, 800_000, 900_000)
        assert isinstance(iterator, abc.Iterable), "returns an iterator"
        all_records = list(iterator)
        assert len(all_records) == 1, "one record in region"

        # Then try a different region. Should reuse same iterator.
        iter2 = simple_tabix_reader.fetch("X", 2_600_000, 3_000_000)
        all_records = list(iter2)
        assert len(all_records) == 3, "Several records in region"

    def test_throws_an_error_if_index_not_present(self):
        with pytest.raises(FileNotFoundError):
            reader = readers.TabixReader(filename=os.path.join(os.path.dirname(__file__), "data/unsorted.txt"),
                                         parser=parsers.standard_gwas_parser)

            reader.fetch('1', 2, 3)

    def test_filter_criteria_with_tabix(self, simple_tabix_reader):
        iterator = simple_tabix_reader.add_filter('ref', 'G').fetch("X", 2_600_000, 3_000_000)
        all_records = list(iterator)
        assert len(all_records) == 2, "The region query obeyed filters"


class TestFileReader:
    def test_filemode_iterator(self, simple_file_reader):
        assert isinstance(simple_file_reader, abc.Iterable)

        iterator = iter(simple_file_reader)
        first_row = next(iterator)
        assert isinstance(first_row, tuple), "rows are parsed as tuples"

    def test_must_specify_filename_or_source(self):
        with pytest.raises(Exception):
            readers.TextFileReader(filename="anything", source=[])

    def test_can_specify_iterable_as_source(self):
        # Any row-based iterable can serve as the source and will be parsed as expected
        parser = readers.TextFileReader(source=["X\t1\tA\tG"])
        result = next(iter(parser))
        assert result[0] == "X"

    def test_can_iterate_twice_over_same_file(self, simple_file_reader):
        iterator = iter(simple_file_reader)
        first_row = next(iterator)
        assert first_row[0] == "1", "Read first data row on first iteration"

        iterator = iter(simple_file_reader)
        first_row = next(iterator)
        assert first_row[0] == "1", "Read first data row on second iteration"


class TestFiltering:
    def test_filter_criteria_with_method(self, simple_file_reader):
        # Can define a filter using any function
        simple_file_reader.add_filter("chrom", lambda val, row: val == "1")
        assert len(simple_file_reader._filters) == 1

        # File will act on it
        assert len(list(simple_file_reader)) == 7, "output was restricted to the expected rows"

    def test_filter_criteria_with_value(self, simple_file_reader):
        # Can define a filter that exactly matches a value
        simple_file_reader.add_filter("chrom", "1")
        assert len(simple_file_reader._filters) == 1

        # File will act on it
        assert len(list(simple_file_reader)) == 7, "output was restricted to the expected rows"

    def test_can_filter_by_tuple_index(self, simple_file_reader):
        # Can apply filters even when the parsed output does not support accessing fields by name
        simple_file_reader.add_filter(0, "1")  # because chr is the first field in the named tuple
        assert len(simple_file_reader._filters) == 1

        # File will act on it
        assert len(list(simple_file_reader)) == 7, "output was restricted to the expected rows"
