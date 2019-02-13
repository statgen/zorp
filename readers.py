"""
Reader objects that handle different types of data
"""
import abc
import collections.abc
import gzip
import logging
import os
import typing

import pysam

import exceptions
import parsers

logger = logging.getLogger(__name__)


class BaseReader(abc.ABC):
    """Implements common base functionality for reading and filtering GWAS results"""
    # TODO: Add a mechanism to skip N rows on iteration, to handle headers vs data
    def __init__(self,
                 filename: str = None,
                 source: typing.Iterable[str] = None,
                 parser: typing.Callable = None,
                 skip_rows: int = 0,
                 skip_errors: bool = False,
                 **kwargs):
        if filename is not None and source is not None:
            # The source option is especially useful for unit testing: create a reader from any iterable
            raise Exception("filename and source options are mutually exclusive")

        # Options for full iteration over file contents (or manually provide a source, which is nice for unit testing)
        self._filename = filename
        self._source = source

        # Filters and transforms that should operate on each row
        # If no parser is provided, just splits the row into a tuple based on the specified delimiter
        self._parser = parser or parsers.TupleLineParser()
        self._filters = []

        self._skip_rows = skip_rows

        # If using "skip error" mode, store a record of which lines had a problem
        self._skip_errors = skip_errors
        self.errors = []

    ######
    # Internal helper methods; should not need to override below this line
    def _apply_filters(self, parsed_row: tuple) -> bool:
        """Determine whether a given row satisfies all of the applied filters"""
        return all(test_func(parsed_row[position], parsed_row)
                   for position, test_func in self._filters)

    @abc.abstractmethod
    def _create_iterator(self) -> typing.Iterable:
        """Create an iterator that can be used over all possible data. Eg, open a file or go through a list"""
        raise NotImplementedError

    def _make_generator(self, iterator) -> typing.Iterable[tuple]:
        """
        Process the output of an iterable before returning it. Omit rows that could not be processed, or that do not
        match the filter criteria.
        """
        for row in iterator:
            try:
                parsed = self._parser(row)
            except exceptions.ParseError as e:
                self.errors.append(row)
                if not self._skip_errors:
                    raise e
                continue

            if not self._apply_filters(parsed):
                continue
            yield parsed

    ######
    # User-facing API
    def add_filter(self,
                   field_name: typing.Union[int, str],
                   match: typing.Union[typing.Any,
                                       typing.Callable[[typing.Any, tuple], bool]]) -> 'BaseReader':
        """
        Limit the output to rows that match the specified criterion. Can apply multiple filters.

        Since the parser should return a namedtuple, can access the field either by name, or by column number

        `match` specified the criterion used to test this field/row. It can be:
         - A value that will exactly match, or
         - A function that specifies whether to accept this row. Method signature: (val, row) => bool
        """
        if not hasattr(self._parser, 'fields'):
            raise Exception('The specified parser does not support accessing properties by name')

        if isinstance(field_name, int):
            position = field_name
        elif field_name not in self._parser.fields:
            raise Exception(f'Parser does not define a field named {field_name}')
        else:
            position = self._parser.fields.index(field_name)

        self._filters.append([
            position,
            match if isinstance(match, collections.abc.Callable) else lambda val, row: val == match
        ])
        return self

    def __iter__(self):
        """
        The parser instance can be treated as an iterable
        """
        # Reset the error list on any new iteration
        self.errors = []

        if self._source:
            return self._make_generator(self._source)
        else:
            iterator = self._create_iterator()
            # Advance the iterator (if applicable) so that only data is returned (not headers)
            for n in range(self._skip_rows):
                next(iterator)
            return self._make_generator(iterator)


class TextFileReader(BaseReader):
    """Implement additional helper methods unique to working with uncompressed text files"""
    def _create_iterator(self):
        """Open the file for parsing"""
        with open(self._filename, 'r') as f:
            yield from f


class TabixReader(BaseReader):
    """
    A reader that can handle generic gzipped data. If a `basefilename.tbi` index is present in the same folder,
    it unlocks additional region-based fetch features.
    """
    def __init__(self, *args, **kwargs):
        self._tabix = None
        super(TabixReader, self).__init__(*args, **kwargs)

        filename = kwargs.get('filename')
        self._has_index = bool(filename and os.path.isfile(f'{filename}.tbi'))

    def _create_iterator(self):
        """
        Return an iterator over all rows of the file
        """
        with gzip.open(self._filename, 'rt') as f:
            yield from f

    def fetch(self, chrom: str, start: int, end: int) -> typing.Iterable:
        """Fetch data from the file in a specific region, by wrapping tabix functionality to return parsed data."""
        if not self._filename:
            raise NotImplementedError("Must specify filename when using region-based fetch")
        if not self._has_index:
            raise FileNotFoundError("You must generate a tabix index before using region-based fetch")

        if not self._tabix:
            # Allow repeatedly fetching from the same tabix file
            self._tabix = pysam.TabixFile(self._filename)

        iterator = self._tabix.fetch(chrom, start, end)
        return self._make_generator(iterator)
