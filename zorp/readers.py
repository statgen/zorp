"""
Reader objects that handle different types of data
"""
import abc
import collections.abc
import gzip
import logging
import os
import typing as ty

import pysam

from . import (
    exceptions,
    parsers,
)


logger = logging.getLogger(__name__)


class BaseReader(abc.ABC):
    """Implements common base functionality for reading and filtering GWAS results"""
    def __init__(self,
                 source: ty.Any,
                 parser: ty.Union[ty.Callable, None] = parsers.TupleLineParser(),
                 skip_rows: int = 0,
                 skip_errors: bool = False,
                 max_errors: int = 100,
                 **kwargs):
        """
        :param source: How to find the source data (implementation detail of the subclass)
        :param parser: A callable that parses a line of text. If `None`, the reader will not attempt to parse.
        :param skip_rows: How many header rows to skip when iterating over the entire file
        :param skip_errors: Whether to continue reading if a few bad lines are encountered
        :param max_errors: How many bad lines to accept (in skip_errors mode) before giving up
        :param kwargs:
        """
        self._source = source

        # Filters and transforms that should operate on each row
        # If no parser is provided, just splits the row into a tuple based on the specified delimiter
        self._parser = parser
        self._filters: list = []

        # If using "skip error" mode, store a record of which lines had a problem (up to a point)
        self._skip_errors = skip_errors
        self._max_errors = max_errors
        self.errors: list = []

        # The user must specify how many header rows to skip
        self._skip_rows = skip_rows

    @abc.abstractmethod
    def _create_iterator(self) -> ty.Iterator[str]:
        """Create an iterator that can be used over all possible data. Eg, open a file or go through a list"""
        raise NotImplementedError

    ######
    # Internal helper methods; should not need to override below this line
    def _apply_filters(self, parsed_row: tuple) -> bool:
        """Determine whether a given row satisfies all of the applied filters"""
        return all(test_func(getattr(parsed_row, field_name, None), parsed_row)
                   for field_name, test_func in self._filters)

    def _make_generator(self, iterator: ty.Iterator[str]) -> ty.Iterator:
        """
        Process the output of an iterable before returning it. In the usual case where a parser is provided,
            this will parse the row, apply filter criteria, and exclude any rows with errors (though errors will be
            available for inspection later)
        """
        use_filters = len(self._filters)
        for i, row in enumerate(iterator):
            if not row:
                # Skip blank lines (eg at end of file)
                continue

            if not self._parser:
                # There is a "parser=None" option, to return raw lines of text. This is useful for, eg, format sniffers.
                yield row
            else:
                try:
                    parsed = self._parser(row)
                except exceptions.LineParseException as e:
                    if not self._skip_errors:
                        raise e
                    self.errors.append((i + self._skip_rows + 1, str(e), row))  # (human_line, message, raw_input)
                    if len(self.errors) >= self._max_errors:
                        raise exceptions.TooManyBadLinesException(error_list=self.errors)
                    continue

                if use_filters and not self._apply_filters(parsed):  # avoid expensive function call if no filters used
                    continue
                yield parsed

    ######
    # User-facing API
    def add_filter(self,
                   field_name: ty.Union[int, str],
                   match: ty.Union[ty.Any,
                                   ty.Callable[[ty.Any, object], bool]]) -> 'BaseReader':
        """
        Limit the output to rows that match the specified criterion. Can apply multiple filters.

        The `field_name` should match a field in whatever object the parser returns. (attribute, calc'd property, etc)

        `match` specified the criterion used to test this field/row. It can be:
         - A value that will exactly match, or
         - A function that specifies whether to accept this row. Method signature: (val, row) => bool
        """
        if not self._parser:
            raise exceptions.ConfigurationException(
                "Filtering features require specifying a parser that supports name-based field access.")

        # Sanity check if possible, otherwise, just hope the parser returns something with named fields!
        if hasattr(self._parser, 'fields') and field_name not in self._parser.fields:  # type: ignore
            raise exceptions.ConfigurationException("The parser does not have a field by this name")

        self._filters.append([
            field_name,
            match if isinstance(match, collections.abc.Callable) else lambda val, row: val == match  # type: ignore
        ])
        return self

    def write(self,
              out_fn: str, *,
              columns: ty.Iterable[str] = None,
              delimiter: str = '\t',
              make_tabix: bool = False):
        """
        Write an output file with the specified columns. Column names must match the field names or tuple indices
        used by your parser.

        This method is useful when you want to save your filtered results for later, or normalize different input
        formats into a single predictable output format.

        FIXME: In tabix mode, out_fn != result_fn . That is... just needlessly confusing.

        TODO: In the initial version, we hardcode a preset mode of "vcf" for tabix indexing.
            There's no good parser-agnostic way to set the tabix options, so we may just have to proxy the pysam kwargs

        # TODO: This isn't a universal converter, eg the column headers are still dictated by the parser (not the user)
        """
        if make_tabix:
            delimiter = '\t'

        if isinstance(self._source, str) and out_fn == self._source:
            raise exceptions.ConfigurationException('Writer cannot overwrite input file')

        if columns is None:
            if hasattr(self._parser, 'fields'):
                columns = self._parser.fields  # type: ignore
            else:
                raise exceptions.ConfigurationException('Must provide column names to write')

        has_headers = any(isinstance(field, str) for field in columns)

        # Special case rule: The writer renders missing data (the Python value `None`) as `.`
        def repr_missing(v):
            return '.' if v is None else v

        with open(out_fn, 'w') as f:
            if has_headers:
                f.write('#')
                f.write(delimiter.join(str(name) for name in columns))
                f.write('\n')
            for row in self:
                f.write(delimiter.join(str(repr_missing(getattr(row, field_name)))
                                       for field_name in columns))
                f.write('\n')

        if make_tabix:
            return pysam.tabix_index(out_fn, force=True, preset='vcf')
        else:
            return out_fn

    def __iter__(self) -> ty.Iterator:
        """
        The parser instance can be treated as an iterable over raw or parsed lines (based on options)
        """
        # Reset the error list on any new iteration
        self.errors = []

        iterator = self._create_iterator()
        # Advance the iterator (if applicable) so that only data is returned (not headers)
        for n in range(self._skip_rows):
            next(iterator)
        return self._make_generator(iterator)


class IterableReader(BaseReader):
    """
    Read data from an externally provided iterable object (like piping from sys.stdin, or a list)

    This can be used as an adapter to wrap an external source in a standard interface, and is also useful for
        unit testing.
    """
    def _create_iterator(self):
        return iter(self._source)


class TextFileReader(BaseReader):
    """Implement additional helper methods unique to working with uncompressed text files"""
    def __init__(self, *args, **kwargs):
        super(TextFileReader, self).__init__(*args, **kwargs)

    def _create_iterator(self):
        """Open the file for parsing"""
        with open(self._source, 'r') as f:
            yield from f


class TabixReader(BaseReader):
    """
    A reader that can handle any gzipped text file. If a `filename.tbi` index is present in the same folder,
    it unlocks additional region-based fetch features.
    """
    def __init__(self, *args, **kwargs):
        super(TabixReader, self).__init__(*args, **kwargs)
        self._tabix: pysam.TabixFile = None
        self._has_index = bool(self._source and os.path.isfile(f'{self._source}.tbi'))

    def _create_iterator(self):
        """
        Return an iterator over all rows of the file
        """
        with gzip.open(self._source, 'rt') as f:
            yield from f

    def fetch(self, chrom: str, start: int, end: int) -> ty.Iterable:
        """Fetch data from the file in a specific region, by wrapping tabix functionality to return parsed data."""
        if not self._has_index:
            raise FileNotFoundError("You must generate a tabix index before using region-based fetch")

        if not self._tabix:
            # Allow repeatedly fetching from the same tabix file
            self._tabix = pysam.TabixFile(self._source)

        iterator = self._tabix.fetch(chrom, start, end)
        return self._make_generator(iterator)


def standard_gwas_reader(filename: str, *, parser=parsers.standard_gwas_parser_quick, **kwargs):  # pragma: no cover
    """
    Helper: generate a reader based on zorp's "standard" format of known column names and positions

    By default this uses the fast (but fragile) "quick" parser, which is not suitable unless the data exactly follows
        the standard format. This may be replaced with a more tolerant parser if needed.
    """
    return TabixReader(filename, parser=parser, skip_rows=1, **kwargs)
