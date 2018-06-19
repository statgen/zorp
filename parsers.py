"""
File format parsers
"""

import abc
import logging
import tempfile
import typing

import pysam


logger = logging.getLogger(__name__)


# TODO: What are some common operations we might want to perform on these files? (alignment across readers, standardization, etc)
# TODO: Standardize the data format (chr, pos, ref, alt, extra
# TODO: Possibly add a native sorting functionality? (eg temp file based?) How to deal with str/unknown chr names?


# TODO: Implement a metaclass that: 1. builds the namedtuple returned type class from fields


all_parsers = []


class BaseGWASMeta(type):
    """
    Metaclass coordinates a registry of known parsers
    """
    def __new__(cls, *args, **kwargs):
        new_class = super(BaseGWASMeta, cls).__new__(cls, *args, **kwargs)
        all_parsers.append(new_class)
        return new_class


class BaseGWASParser(metaclass=BaseGWASMeta):
    def __init__(self,
                 filename: str=None,
                 source: typing.Iterable[str]=None,
                 has_headers: bool=True,
                 **kwargs):
        if filename is not None and source is not None:
            raise Exception("filename and source options are mutually exclusive")

        # Options for full iteration over file contents
        self.filename = filename
        self._iterator = source

        # Options for tabix region iteration
        self._tabix = None

        # Filters and transforms
        self._filters = []
        self._transforms = []  # TODO: Implement transforms

    ######
    # Custom behavior that must be implemented in subclasses
    def _parse_row(self, row: str) -> tuple:
        """Handle the logic of parsing data from the specific file format"""
        # TODO: Which columns? Do values require any cleanup, like chr/pos/ref/alt?
        raise NotImplementedError
        # return tuple()


    ######
    # User-facing API
    def fetch(self, chr:str, start: int, end: int) -> typing.Iterable:
        """Fetch data from the file in a specific region, by wrapping tabix functionality to return parsed data."""
        if not self.filename:
            raise Exception("Must specify filename when using region-based fetch")

        if not self._tabix:
            # Allow repeatedly fetching from the same tabix file
            self._tabix = pysam.TabixFile(self.filename)

        return self._make_generator(self._tabix.fetch(chr, start, end))

    def add_filter(self,
                   field_name: str,
                   test_func: typing.Callable[[typing.Any, tuple], bool]) -> 'BaseGWASParser':
        """
        Limit the output to rows that match the specified criterion. Can apply multiple filters.
        """
        # TODO: Validate that we are filtering a known field?
        self._filters.append([field_name, test_func])
        return self

    def sort(self, lazy: bool=True) -> 'BaseGWASParser':
        """
        Sort the file provided, and return a new instance with that handle as the iterator
        :param lazy: For large files, this will use a simple heuristic (smaller of first 10% or thousand lines)
            to decide whether sorting is necessary.
        """

        if self._iterator:
            raise Exception("Sorting only available on files")

        # TODO: probably an instance method, not a class... and should check file isn't open yet, set its own iterator to file handle
        with tempfile.TemporaryFile() as f:  # TODO: Sort, write, seek to 0
            pass
        f.seek(0)

        raise NotImplementedError("Must implement sorting logic for this file type")
        return self

    ######
    # Internal helper methods; should not need to override below this line
    def _apply_filters(self, parsed_row: tuple) -> bool:
        """Determine whether a given row satisfies all of the applied filters"""
        return all(test_func(getattr(parsed_row, field_name), parsed_row)
                   for field_name, test_func in self._filters)

    def _create_file_iterator(self):
        """Open the file for parsing"""
        with open(self.filename) as f:
            yield from f

    def _make_generator(self, iterator) -> typing.Iterable:
        """Process the output of an iterable before returning it"""
        for row in iterator:
            parsed = self._parse_row(row)
            if not self._apply_filters(parsed):
                continue
            yield parsed

    def __iter__(self):
        """
        The parser instance can be treated as an iterable
        """
        if not self._iterator:
            self._iterator = self._create_file_iterator()

        return self._make_generator(self._iterator)

# Write parsers for raremetal,rvtests, possibly others. See edge cases peter handles and crib mercilessly


