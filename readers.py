"""
File format parsers
"""
import abc
import collections.abc
import gzip
import logging
import subprocess
import tempfile
import typing

import pysam


logger = logging.getLogger(__name__)

# TODO: What are some common operations we might want to perform on these files? (alignment across readers, standardization, etc)
#   Ensure that base classes are provided to support this


class BaseReader(abc.ABC):
    """Implements common base functionality for reading and filtering GWAS results"""
    def __init__(self,
                 filename: str=None,
                 source: typing.Iterable[str]=None,
                 parser: typing.Callable=None,
                 **kwargs):
        if filename is not None and source is not None:
            raise Exception("filename and source options are mutually exclusive")

        # Options for full iteration over file contents
        self.filename = filename
        self._iterator = source

        # Filters and transforms that should operate on each row
        # If no parser is provided, just splits the row into a tuple based on the specified delimiter
        self._parser = parser or (lambda row: tuple(row.split(kwargs.get('delimiter', '\t'))))
        self._filters = []

    ######
    # User-facing API
    def add_filter(self,
                   field_name: str,
                   match: typing.Union[typing.Any,
                                       typing.Callable[[typing.Any, tuple], bool]]) -> 'BaseReader':
        """
        Limit the output to rows that match the specified criterion. Can apply multiple filters.

        `match` specified the criterion used to test this field/row. It can be:
         - A value that will exactly match, or
         - A function that specifies whether to accept this row. Method signature: (val, row) => bool
        """
        if not hasattr(self._parser, 'fields'):
            raise Exception('The specified parser does not support accessing properties by name')

        if field_name not in self._parser.fields:
            raise Exception(f'Parser does not define a field named ${field_name}')

        self._filters.append([
            field_name,
            match if isinstance(match, collections.abc.Callable) else lambda val, row: val == match
        ])
        return self

    ######
    # Internal helper methods; should not need to override below this line
    def _apply_filters(self, parsed_row: tuple) -> bool:
        """Determine whether a given row satisfies all of the applied filters"""
        return all(test_func(getattr(parsed_row, field_name), parsed_row)
                   for field_name, test_func in self._filters)

    @abc.abstractmethod
    def _create_iterator(self) -> typing.Iterable:
        """Create an iterator that can be used over all possible data. Eg, open a file or go through a list"""
        raise NotImplementedError

    def _make_generator(self, iterator) -> typing.Iterable:
        """Process the output of an iterable before returning it"""
        for row in iterator:
            parsed = self._parser(row)
            if not self._apply_filters(parsed):
                continue
            yield parsed

    def __iter__(self):
        """
        The parser instance can be treated as an iterable
        """
        if not self._iterator:
            self._iterator = self._create_iterator()

        return self._make_generator(self._iterator)


class FileReader(BaseReader):
    """Implement additional helper methods unique to working with uncompressed text files"""
    def sort(self, primary_col, secondary_col, skip_rows=0) -> 'FileReader':
        """
        Sort the file provided, and update the internal iterator
        """
        if self._iterator:
            raise Exception("Sorting only available on unopened files")

        # TODO: handle skipping headers when necessary
        f = tempfile.TemporaryFile(mode='r+')
        subprocess.check_call([
            'sort',
            self.filename,
            '-k{0},{0}'.format(primary_col),  # chr can be string (or not; we want some numeric sorting)
            '-k{0},{0}n'.format(secondary_col)  # pos should be numeric
        ], stdout=f)
        f.seek(0)
        self._iterator = f
        return self

    def _create_iterator(self):
        """Open the file for parsing"""
        with open(self.filename) as f:
            yield from f


class TabixReader(BaseReader):
    def __init__(self, *args, **kwargs):
        self._tabix = None
        super(TabixReader, self).__init__(*args, **kwargs)

    def _create_iterator(self):
        """
        Return an iterator over all rows of the file
        """
        with gzip.open(self.filename, 'rt') as f:
            yield from f

    def fetch(self, chr:str, start: int, end: int) -> typing.Iterable:
        """Fetch data from the file in a specific region, by wrapping tabix functionality to return parsed data."""
        if not self.filename:
            raise Exception("Must specify filename when using region-based fetch")

        if not self._tabix:
            # Allow repeatedly fetching from the same tabix file
            self._tabix = pysam.TabixFile(self.filename)

        self._iterator = self._tabix.fetch(chr, start, end)
        return self
