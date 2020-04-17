"""
Reader objects that handle different types of data
"""
import abc
import collections.abc
import gzip
import logging
import os
import sys
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
                 parser: ty.Union[ty.Callable[[str], object], None] = parsers.TupleLineParser(),
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

        # Filters, lookups, and mutations that should operate on each row
        # If no parser is provided, just splits the row into a tuple based on the specified delimiter
        self._parser = parser
        self._filters = []  # type: list  # Should we return this row from iteration?
        self._lookups = []  # type: list  # Find the value of a specified field given (parsed) variant info
        self._transforms = []  # type: list  # Modify the (parsed) variant info using a custom function.

        # If using "skip error" mode, store a record of which lines had a problem (up to a point)
        self._skip_errors = skip_errors
        self._max_errors = max_errors
        self.errors = []  # type: list

        # The user must specify how many header rows to skip
        self._skip_rows = skip_rows

    @abc.abstractmethod
    def _create_iterator(self) -> ty.Iterator[str]:
        """Create an iterator that can be used over all possible data. Eg, open a file or go through a list"""
        raise NotImplementedError

    ######
    # Internal helper methods; should not need to override below this line

    def _make_generator(self, iterator: ty.Iterator[str]) -> ty.Iterator:
        """
        Process the output of an iterable before returning it. In the usual case where a parser is provided,
            this will parse the row, apply filter criteria, and exclude any rows with errors (though errors will be
            available for inspection later)
        """
        use_filters = len(self._filters)
        use_lookups = len(self._lookups)
        use_transforms = len(self._transforms)
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

                if use_lookups:
                    for field_name, func in self._lookups:
                        setattr(parsed, field_name, func(parsed))

                if use_transforms:
                    for func in self._transforms:
                        parsed = func(parsed)

                if use_filters and not all(test_func(parsed) for test_func in self._filters):
                    continue
                yield parsed

    ######
    # User-facing API
    def add_filter(self, *args) -> 'BaseReader':
        """
        Limit the output to rows that match the specified criterion. Can apply multiple filters.
        Filters are applied after all parsing and transforms are complete.

        There are three ways to specify a filter:
        - `add_filter('field_name') requires that the value not be missing. Equivalent to "is not None".
        - `add_filter('field_name', value)` checks that the field has this exact value.
        - `add_filter(lambda parsed: bool)` runs a user-provided function on whatever is in this row.

        The `field_name` should match a field in whatever object the parser returns. (attribute, calc'd property, etc)
        """
        if len(args) == 1:
            spec = args[0]
            if isinstance(spec, collections.abc.Callable):  # type: ignore
                self._filters.append(spec)
            elif isinstance(spec, str):
                self._filters.append(lambda parsed: getattr(parsed, spec) is not None)
            else:
                raise exceptions.ConfigurationException('Single argument must be either a function or a field name')
        elif len(args) == 2:
            # Exact value match
            field_name, target_value = args
            self._filters.append(lambda parsed: getattr(parsed, field_name) == target_value)
        else:
            raise exceptions.ConfigurationException('Invalid filter format requested')
        return self

    def add_lookup(self, field_name: str, lookup_func: ty.Callable[[object], object]) -> 'BaseReader':
        """
        Look up / modify the value of an individual field within a row. Each lookup is a function that receives
          the parsed row and returns a single value. This can be used to clean up how values are presented
          ("chr1" vs "chr1"), and also for things that a generic parser may be ill-suited to infer from the text of the
          file. (eg, rsid lookups require knowing reference build, which is usually specified as metadata).

        Lookups are good at finding the value of a single field. Lookups can be written that have no dependence on
          zorp, eg "genelocator".

        Example: `add_lookup('rsid', lookup_func)` sets `parsed.rsid` to the result of the lookup.

        At present, zorp is written with the restrictive philosophy that each field has a well-understood universal
          meaning, consistent across all files that use that parser. Lookups are allowed to modify existing known
          fields, but not add new fields. (a side effect: the parser must support name-based field access)
        """
        if not self._parser:
            raise exceptions.ConfigurationException(
                "Lookup features require specifying a parser that supports name-based field access.")

        # Sanity check if possible, otherwise, just hope the parser returns something with named fields!
        if hasattr(self._parser, 'fields') and field_name not in self._parser.fields:  # type: ignore
            raise exceptions.ConfigurationException("The parser does not have a field by this name")

        if not isinstance(lookup_func, collections.abc.Callable):  # type: ignore
            raise exceptions.ConfigurationException(
                "Lookup must specify a function that receives variant data and returns a value")

        self._lookups.append([
            field_name,
            lookup_func
        ])
        return self

    def add_transform(self, transform_func: ty.Callable[[object], object]) -> 'BaseReader':
        """
        Transforms provide a way to change ALL of the data for a variant. They can mutate several fields at once,
          perform multiple lookups, or even return a different object than what the parser normally returns. (this last
          option is not recommended for both performance and debugging reasons)

        We generally recommend lookups for most use cases, because they don't need to depend on your parser,
          and are easier to maintain. (because they are pure functions with no side effects). Transforms represent
          a way to make very powerful arbitrary changes to all your data at once, but they also bypass type checks and
          other features designed to prevent bugs.
        """
        if not isinstance(transform_func, collections.abc.Callable):  # type: ignore
            raise exceptions.ConfigurationException(
                "Transforms must specify a function that operates on variant data")
        self._transforms.append(transform_func)
        return self

    def write(self,
              out_fn: str = None, *,
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

        if out_fn is None and make_tabix:
            raise exceptions.ConfigurationException('Cannot create tabix file if writing output to a stream')

        if isinstance(self._source, str) and out_fn == self._source:
            raise exceptions.ConfigurationException('Writer cannot overwrite input file')

        if columns is None:
            try:
                columns = self._parser.fields  # type: ignore
            except AttributeError:
                raise exceptions.ConfigurationException('Must provide column names to write')

        # Special case rule: The writer renders missing data (the Python value `None`) as `.`
        def repr_missing(v):
            return '.' if v is None else v

        def write_all(handle):
            """Internal helper that allows writing to either a file, or stdout"""
            # Write headers
            try:
                handle.write('#')
                handle.write(delimiter.join(str(name) for name in columns))
                handle.write('\n')
                for row in self:
                    handle.write(delimiter.join(str(repr_missing(getattr(row, field_name)))
                                                for field_name in columns))
                    handle.write('\n')
            except BrokenPipeError:  # pragma: no cover
                # When writing to stdout, some utils (like head) may close the pipe early, at which point we end writing
                return

        # Readers can write to stdout, which lets CLI scripts (like zorp-convert) use this in a pipeline
        try:
            with open(out_fn, 'w') as f:
                write_all(f)
        except TypeError:
            write_all(sys.stdout)

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
        self._tabix = None  # type: pysam.TabixFile
        self._has_index = bool(self._source and os.path.isfile('{}.tbi'.format(self._source)))

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
