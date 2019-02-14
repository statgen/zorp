"""
Parsers: handle the act of reading one entity (such as line)
"""
import abc
import numbers
import re
import typing

from . import exceptions

# TODO: Improve marker pattern
MARKER_PATTERN = re.compile(r'([^:]+):([0-9]+)_([-ATCG.]+)/([-ATCG.*]+)')


class _basic_standard_container(typing.NamedTuple):
    """Store GWAS results in a predictable format, with only the minimum fields"""
    chrom: str
    pos: int
    ref: str
    alt: str
    pvalue: numbers.Number


class _extended_standard_container(typing.NamedTuple):
    """Store GWAS results in a predictable format, with a variety of added fields"""
    # Required fields
    chrom: str
    pos: int
    ref: str
    alt: str
    pvalue: numbers.Number

    # Optional fields
    af: float
    beta: float
    stderr: float
    marker: str
    rsid: str


class AbstractLineParser(abc.ABC):
    """
    Abstract parser that returns a line of text into GWAS data
    This base class is reserved for future refactoring
    """
    def __init__(self, *args,
                 container: typing.Callable[..., typing.Union[tuple, typing.NamedTuple]] = tuple,
                 **kwargs):
        """
        :param container: A data structure (eg namedtuple) that will be populated with the parsed results
        """
        self._container = container

    @abc.abstractmethod
    def _split_fields(self, row: str):
        """Turn a row of text into separate fields, eg by splitting a delimiter or deserializing JSON"""
        pass

    @abc.abstractmethod
    def _process_values(self, values: typing.Sequence) -> tuple:
        """Convert the raw field values into something useful. This includes field selection and type coercion."""
        pass

    @abc.abstractmethod
    def _output_container(self, values: typing.Iterable):
        """Populate the output container with the extracted results"""
        pass

    def __call__(self, row: str) -> tuple:
        fields = self._split_fields(row)
        values = self._process_values(fields)
        container = self._output_container(values)
        return container


class TupleLineParser(AbstractLineParser):
    """Parse a line of text and return a tuple of the fields. Performs no type coercion"""
    def __init__(self, *args, container: typing.Callable = tuple, delimiter='\t', **kwargs):
        super(TupleLineParser, self).__init__(*args, container=container, **kwargs)
        self._delimiter = delimiter

    def _split_fields(self, row: str):
        return row.strip().split(self._delimiter)

    def _process_values(self, values: typing.Sequence):
        return values

    def _output_container(self, values):
        return self._container(values)


class GenericGwasLineParser(TupleLineParser):
    """
    A simple parser that extracts GWAS information from a flexible file format.

    Column numbers in the parser are 0-delimited
    """
    def __init__(self, *args,
                 container: typing.Callable[..., tuple] = _basic_standard_container,
                 chr_col: int = None, pos_col: int = None, ref_col: int = None, alt_col: int = None,
                 marker_col: int = None, rsid_col: int = None,
                 pval_col: int = None,
                 # Optional fields
                 af_col: int = None, beta_col: int = None, stderr_col: int = None,
                 # Other configuration
                 is_log_pval: bool = False,
                 **kwargs):
        super(GenericGwasLineParser, self).__init__(*args, **kwargs)

        # All GWAS parsers can specify
        self._container = container

        # The kwargs use human-friendly numbers. Internally, we store them as 0-based indices.
        def _human_to_zero(value):
            return value - 1 if value else None

        self._chr_col = _human_to_zero(chr_col)
        self._pos_col = _human_to_zero(pos_col)
        self._ref_col = _human_to_zero(ref_col)
        self._alt_col = _human_to_zero(alt_col)

        self._marker_col = _human_to_zero(marker_col)
        self._rsid_col = _human_to_zero(rsid_col)

        self._pval_col = _human_to_zero(pval_col)

        self._af_col = _human_to_zero(af_col)
        self._beta_col = _human_to_zero(beta_col)
        self._stderr_col = _human_to_zero(stderr_col)

        self._is_log_pval = is_log_pval

    # def _parse_marker(self, marker) -> tuple:
    #     match = MARKER_PATTERN.match(marker)
    #     if match is None:
    #         raise exceptions.ParseError(f'Could not read interpret marker ID: {marker}')
    #
    #     return match.groups()

    @property
    def fields(self) -> typing.Container:
        return self._container._fields  # type: ignore

    def _get_pval(self, pvalue):
        pvalue = float(pvalue)
        if self._is_log_pval:
            pvalue = 10 ** -pvalue
        return pvalue

    def _process_values(self, values: typing.Sequence):
        # Fetch values
        chrom = values[self._chr_col]
        pos = values[self._pos_col]
        ref = values[self._ref_col]
        alt = values[self._alt_col]

        pval = values[self._pval_col]

        # Perform type coercion
        try:
            pval = self._get_pval(pval)
            pos = int(pos)
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=values)

        return chrom, pos, ref, alt, pval

    def _output_container(self, values):
        return self._container(*values)


# An example parser pre-configured for the LocusZoom standard file format
standard_gwas_parser = GenericGwasLineParser(chr_col=1, pos_col=2, ref_col=3, alt_col=4,
                                             pval_col=5, is_log_pval=True,
                                             delimiter='\t')
