"""
Parsers: handle the act of reading one entity (such as line)
"""
import abc
import inspect
import math
import numbers
import typing as ty

from .const import MISSING_VALUES
from . import exceptions, parser_utils


class _basic_standard_container(ty.NamedTuple):
    """Store GWAS results in a predictable format, with only the minimum fields"""
    chrom: str
    pos: int
    ref: str
    alt: str
    neg_log_pvalue: float

    # # Optional fields for future expansion
    # af: float
    # beta: float
    # stderr: float
    # marker: str
    # rsid: str

    @property
    def pvalue(self) -> ty.Union[float, None]:
        if self.neg_log_pvalue is None:
            return None
        elif math.isinf(self.neg_log_pvalue):
            # This is an explicit design choice here, since we parse p=0 to infinity
            return 0
        else:
            return 10 ** -self.neg_log_pvalue

    @property
    def pval(self) -> numbers.Number:
        """A common field name alias"""
        return self.pvalue

    @property
    def marker(self) -> str:
        """Specify the marker in a string format compatible with UM LD server and other variant-specific requests"""
        ref_alt = f'_{self.ref}/{self.alt}' if (self.ref and self.alt) else ''
        return f'{self.chrom}:{self.pos}{ref_alt}'

    @classmethod
    def _to_serialize(cls)-> ty.List[str]:
        return [name for name, value in inspect.getmembers(cls) if inspect.isdatadescriptor(value)]

    def to_dict(self):
        # A special version of asdict that also includes derived properties
        return {name: getattr(self, name) for name in self._to_serialize()}


class AbstractLineParser(abc.ABC):
    """
    Abstract parser that returns a line of text into GWAS data
    This base class is reserved for future refactoring
    """
    def __init__(self, *args,
                 container: ty.Callable[..., ty.Union[tuple, ty.NamedTuple]] = tuple,
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
    def _process_values(self, values: ty.Sequence) -> tuple:
        """Convert the raw field values into something useful. This includes field selection and type coercion."""
        pass

    @abc.abstractmethod
    def _output_container(self, values: ty.Iterable):
        """Populate the output container with the extracted results"""
        pass

    def __call__(self, row: str) -> tuple:
        try:
            fields = self._split_fields(row)
            values = self._process_values(fields)
            container = self._output_container(values)
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=row)
        return container


class TupleLineParser(AbstractLineParser):
    """Parse a line of text and return a tuple of the fields. Performs no type coercion"""
    def __init__(self, *args, container: ty.Callable = tuple, delimiter='\t', **kwargs):
        super(TupleLineParser, self).__init__(*args, container=container, **kwargs)
        self._delimiter = delimiter

    def _split_fields(self, row: str):
        return row.strip().split(self._delimiter)

    def _process_values(self, values: ty.Sequence):
        return values

    def _output_container(self, values):
        return self._container(values)


class GenericGwasLineParser(TupleLineParser):
    """
    A simple parser that extracts GWAS information from a flexible file format.

    Constructor expects human-friendly column numbers (first = column 1)
    """
    def __init__(self, *args,
                 container: ty.Callable[..., tuple] = _basic_standard_container,
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

        self.validate_config()

    def validate_config(self):
        """Ensures that a minimally working parser has been created"""
        has_position = (self._marker_col is not None) ^ all(getattr(self, x) is not None
                                                            for x in ('_chr_col', '_pos_col', '_ref_col', '_alt_col'))
        is_valid = has_position and (self._pval_col is not None)
        if not is_valid:
            raise exceptions.ConfigurationException('GWAS parser must specify how to find all required fields')
        return is_valid

    @property
    def fields(self) -> ty.Container:
        return self._container._fields  # type: ignore

    def _process_values(self, values: ty.Sequence):
        # Fetch values
        if self._marker_col is not None:
            chrom, pos, ref, alt = parser_utils.parse_marker(values[self._marker_col])
        else:
            # TODO: Should we check for, and strip, the letters chr?
            chrom = values[self._chr_col]
            pos = values[self._pos_col]
            ref = values[self._ref_col]
            alt = values[self._alt_col]

        pval = values[self._pval_col]

        # Perform type coercion
        try:
            log_pval = parser_utils.parse_pval_to_log(pval, is_log=self._is_log_pval)
            pos = int(pos)
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=values)

        # Some old GWAS files simply won't provide ref or alt information, and the parser will need to do without
        if ref in MISSING_VALUES:
            ref = None

        if alt in MISSING_VALUES:
            alt = None

        return chrom, pos, ref, alt, log_pval

    def _output_container(self, values):
        return self._container(*values)


# An example parser pre-configured for the LocusZoom standard file format
standard_gwas_parser = GenericGwasLineParser(chr_col=1, pos_col=2, ref_col=3, alt_col=4,
                                             pval_col=5, is_log_pval=True,
                                             delimiter='\t')
