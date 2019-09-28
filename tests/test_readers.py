"""
Test reader functionality
"""
from collections import abc
import os

import pytest

from zorp import (
    exceptions,
    readers,
)


@pytest.fixture
def simple_tabix_reader(standard_gwas_parser_basic):
    return readers.TabixReader(os.path.join(os.path.dirname(__file__), "data/sample.gz"),
                               parser=standard_gwas_parser_basic,
                               skip_rows=1)


@pytest.fixture
def simple_file_reader(standard_gwas_parser_basic):
    return readers.TextFileReader(os.path.join(os.path.dirname(__file__),
                                               "data/pheweb-samples/has-fields-.txt"),
                                  parser=standard_gwas_parser_basic,
                                  skip_rows=1)


def doomed_parser(line):
    """Fixture to test error handling"""
    raise exceptions.LineParseException('Error occurred')


class TestIterableReader:
    def test_can_specify_iterable_as_source(self):
        # Any row-based iterable can serve as the source and will be parsed as expected
        reader = readers.IterableReader(["X\t1\tA\tG"])
        result = next(iter(reader))
        assert result[0] == "X"

    #####
    # Special cases: Things we should not try to parse
    def test_skips_empty_rows_padding_file(self):
        reader = readers.IterableReader(["", ""])
        results = list(reader)
        assert len(results) == 0, "Skipped empty lines"

    def test_can_optionally_iterate_sans_parsing(self):
        reader = readers.IterableReader(["walrus", "carpenter"], parser=None)
        results = list(reader)
        assert results == ["walrus", "carpenter"], "Returns unparsed data"

    ######
    # Writes files
    def test_can_write_output(self, tmpdir, standard_gwas_parser_basic):
        reader = readers.IterableReader(["1\t100\tA\tC\t0.05", "2\t200\tA\tC\t5e-8"],
                                        parser=standard_gwas_parser_basic)
        expected_fn = tmpdir / 'test.txt'
        out_fn = reader.write(expected_fn, columns=['chrom'], make_tabix=False)

        assert expected_fn == out_fn
        assert os.path.isfile(out_fn), "Output filename exists"
        with open(out_fn, 'r') as f:
            assert f.readlines() == ["#chrom\n", "1\n", "2\n"]

    def test_writer_represents_missing_data_correctly(self, tmpdir, standard_gwas_parser_basic):
        """The writer should represent explicit missing values as `.` (instead of eg Python None)"""
        reader = readers.IterableReader(["1\t100\tA\tC\tNone", "2\t200\tA\tC\t."],
                                        parser=standard_gwas_parser_basic)
        expected_fn = tmpdir / 'test.txt'
        out_fn = reader.write(expected_fn, columns=['neg_log_pvalue'], make_tabix=False)

        assert expected_fn == out_fn
        assert os.path.isfile(out_fn), "Output filename exists"
        with open(out_fn, 'r') as f:
            assert f.readlines() == ["#neg_log_pvalue\n", ".\n", ".\n"]

    def test_can_write_tabixed_output(self, tmpdir, standard_gwas_parser_basic):
        reader = readers.IterableReader(["1\t100\tA\tC\t0.05", "2\t200\tA\tC\t5e-8"],
                                        parser=standard_gwas_parser_basic)
        expected_fn = tmpdir / 'test.gz'
        out_fn = reader.write(str(expected_fn), columns=['chrom', 'pos'], make_tabix=True)

        assert expected_fn != out_fn
        assert out_fn.endswith('.gz')
        assert os.path.exists(f'{out_fn}.tbi'), "Tabix index exists"

        assert os.path.isfile(out_fn), "Output filename exists"

        # Now try to use the file that was written
        check_output = readers.TabixReader(out_fn)
        assert len(list(check_output.fetch('1', 1, 300))) == 1, 'Output file can be read with tabix features'

    def test_writer_defaults_to_parser_columns(self, tmpdir, standard_gwas_parser_basic):
        reader = readers.IterableReader(['1\t100\tA\tC\t0.05', '2\t200\tA\tC\t5e-8'],
                                        parser=standard_gwas_parser_basic)
        expected_fn = tmpdir / 'test.txt'
        out_fn = reader.write(expected_fn)

        with open(out_fn, 'r') as f:
            assert f.readline() == '#chrom\tpos\trsid\tref\talt\tneg_log_pvalue\tbeta\tstderr_beta\talt_allele_freq\n'

    def test_writer_needs_to_know_column_names(self, tmpdir):
        reader = readers.IterableReader(['1\t100\tA\tC\t0.05', '2\t200\tA\tC\t5e-8'],
                                        parser=lambda line: ('A', 'B'))
        expected_fn = tmpdir / 'test.txt'
        with pytest.raises(exceptions.ConfigurationException, match='column names'):
            reader.write(expected_fn)

    def test_writer_can_send_to_console_stdout(self, capsys, standard_gwas_parser_basic):
        reader = readers.IterableReader(['1\t100\tA\tC\t0.05', '2\t200\tA\tC\t5e-8'],
                                        parser=standard_gwas_parser_basic)
        reader.write()
        out, err = capsys.readouterr()
        assert out.splitlines()[0] == '#chrom\tpos\trsid\tref\talt\tneg_log_pvalue\tbeta\tstderr_beta\talt_allele_freq'

    def test_writer_validates_options_when_sending_to_console(self, standard_gwas_parser_basic):
        reader = readers.IterableReader(['1\t100\tA\tC\t0.05', '2\t200\tA\tC\t5e-8'],
                                        parser=standard_gwas_parser_basic)
        with pytest.raises(exceptions.ConfigurationException, match='stream'):
            reader.write(make_tabix=True)

    ######
    # Error handling
    def test_can_fail_on_first_error(self):
        reader = readers.IterableReader(['mwa', 'ha', 'ha'], parser=doomed_parser, skip_errors=False)
        with pytest.raises(exceptions.LineParseException):
            list(reader)

    def test_can_track_errors(self):
        reader = readers.IterableReader(['mwa', 'ha', 'ha'], parser=doomed_parser, skip_errors=True, max_errors=10)
        results = list(reader)
        assert len(results) == 0, "No data could actually be read!"
        assert len(reader.errors) == 3, "Three lines could not be parsed"

    def test_warns_if_file_is_unreadable(self):
        reader = readers.IterableReader(['mwa', 'ha', 'ha'], parser=doomed_parser, skip_errors=True, max_errors=2)
        with pytest.raises(exceptions.TooManyBadLinesException):
            list(reader)
        assert len(reader.errors) == 2, "Reader gave up after two lines, but tracked the errors"


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

    def test_throws_an_error_if_index_not_present(self, standard_gwas_parser):
        with pytest.raises(FileNotFoundError):
            reader = readers.TabixReader(os.path.join(os.path.dirname(__file__), "data/unsorted.txt"),
                                         parser=standard_gwas_parser)
            reader.fetch('1', 2, 3)

    def test_filter_criteria_with_tabix(self, simple_tabix_reader):
        iterator = simple_tabix_reader\
            .add_filter('ref', 'G')\
            .fetch("X", 2_600_000, 3_000_000)
        all_records = list(iterator)
        assert len(all_records) == 2, "The region query obeyed filters"

    def test_can_iterate_over_all_rows(self, simple_tabix_reader):
        assert len(list(simple_tabix_reader)) == 88, 'Fetched all rows from sample file'


class TestFileReader:
    def test_filemode_iterator(self, simple_file_reader):
        assert isinstance(simple_file_reader, abc.Iterable)

        iterator = iter(simple_file_reader)
        first_row = next(iterator)
        assert isinstance(first_row, object), "rows are parsed as tuples"

    def test_can_iterate_twice_over_same_file(self, simple_file_reader):
        iterator = iter(simple_file_reader)
        first_row = next(iterator)
        assert first_row.chrom == "1", "Read first data row on first iteration"

        iterator = iter(simple_file_reader)
        first_row = next(iterator)
        assert first_row.chrom == "1", "Read first data row on second iteration"

    def test_writer_protects_from_overwriting(self, simple_file_reader):
        with pytest.raises(exceptions.ConfigurationException):
            simple_file_reader.write(simple_file_reader._source, columns=[])


class TestFiltering:
    def test_filter_criteria_with_method(self, simple_file_reader):
        # Can define a filter using any function
        simple_file_reader.add_filter(lambda row: row.chrom == "1")
        assert len(simple_file_reader._filters) == 1

        # File will act on it
        assert len(list(simple_file_reader)) == 7, "output was restricted to the expected rows"

    def test_filter_criteria_with_value(self, simple_file_reader):
        # Can define a filter that exactly matches a value
        simple_file_reader.add_filter("chrom", "1")
        assert len(simple_file_reader._filters) == 1

        # File will act on it
        assert len(list(simple_file_reader)) == 7, "output was restricted to the expected rows"

    def test_filter_with_just_field_name_excludes_missing_rows(self, standard_gwas_parser_basic):
        reader = readers.IterableReader(["1\t100\tnull\tNone\t0.05", "2\t200\tA\tC\t5e-8"],
                                        parser=standard_gwas_parser_basic)
        reader.add_filter('ref')
        assert len(reader._filters) == 1

        # File will act on it
        assert len(list(reader)) == 1, "output was restricted to the expected rows"


class TestTransforms:
    def test_transform_modifies_parsed_field_value(self, simple_file_reader):
        reader = simple_file_reader.add_transform('chrom', lambda parsed: 'Y')
        assert all(row.chrom == 'Y' for row in reader), 'All chromosomes have been replaced with a new value'

    def test_named_field_parsers_can_check_that_dest_field_exists(self, simple_file_reader):
        with pytest.raises(exceptions.ConfigurationException, match='does not have a field'):
            simple_file_reader.add_transform('dummy_field', lambda parsed: 'Y')

    def test_transform_must_be_function(self, simple_file_reader):
        with pytest.raises(exceptions.ConfigurationException, match='must specify a function'):
            simple_file_reader.add_transform('chrom', 12)

    def test_named_fields_require_named_field_parser(self):
        reader = readers.IterableReader([[1, 2, 3]], parser=None)
        with pytest.raises(exceptions.ConfigurationException, match='name-based'):
            reader.add_transform('chrom', lambda parsed: 'Y')
