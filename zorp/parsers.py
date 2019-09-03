"""
Parsers: handle the act of reading one entity (such as line)
"""
import abc
import math
import numbers
import typing as ty

try:
    from fastnumbers import int, float
except ImportError:
    pass

from .const import MISSING_VALUES
from . import exceptions, parser_utils


class BasicVariant:
    """
    Store GWAS results in a predictable format, with a minimal set of fields; optimize for name-based attribute access
    """
    # TODO: DRY these lists; this is a bit silly
    # Slots specify the data  this holds (a performance optimization); _fields is human-curated list
    __slots__ = ('chrom', 'pos', 'ref', 'alt', 'neg_log_pvalue', 'beta', 'stderr_beta', 'alt_allele_freq')
    _fields = ('chrom', 'pos', 'ref', 'alt', 'neg_log_pvalue', 'beta', 'stderr_beta', 'alt_allele_freq')

    def __init__(self, chrom, pos, ref, alt, neg_log_pvalue, beta, stderr_beta, alt_allele_freq):
        self.chrom: str = chrom
        self.pos: int = pos
        self.ref: str = ref
        self.alt: str = alt
        self.neg_log_pvalue: float = neg_log_pvalue

        # # Optional fields for future expansion
        # af: float
        self.beta: float = beta
        self.stderr_beta: float = stderr_beta

        self.alt_allele_freq = alt_allele_freq
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

    def __iter__(self):
        # Features like "write a line of text" may want to access fields in a predictable order
        return self.__slots__

    def to_dict(self):
        # Some tools expect the data in a mutable form (eg dicts)
        return {s: getattr(self, s, None) for s in self._fields}


class AbstractLineParser(abc.ABC):
    """
    Abstract parser that returns a line of text into GWAS data
    This base class is reserved for future refactoring
    """
    def __init__(self, *args,
                 container: ty.Callable[..., object] = tuple,
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
    """
    Parse a line of text and return a tuple of the fields. Performs no type coercion

    This isn't recommended for everyday parsing, but it is useful internally for format detection (where we need to
        split columns of data, but aren't yet ready to clean and assign meaning to the values)
    """
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
                 container: ty.Callable[..., BasicVariant] = BasicVariant,
                 chr_col: int = None, pos_col: int = None, ref_col: int = None, alt_col: int = None,
                 marker_col: int = None, rsid_col: int = None,
                 pval_col: int = None,
                 # Optional fields
                 beta_col: int = None, stderr_col: int = None,
                 allele_freq_col: int = None, allele_count_col: int = None, n_samples_col: int = None,
                 # Other configuration options that apply to every row as constants
                 is_log_pval: bool = False,
                 is_alt_effect: bool = True,  # whether effect allele is oriented towards alt
                 **kwargs):
        super(GenericGwasLineParser, self).__init__(*args, **kwargs)

        # All GWAS parsers can specify
        self._container = container

        # The kwargs use human-friendly numbers. Internally, we store them as 0-based indices.
        def _human_to_zero(value):
            if value is None:
                return value
            else:
                return value - 1

        self._chr_col = _human_to_zero(chr_col)
        self._pos_col = _human_to_zero(pos_col)
        self._ref_col = _human_to_zero(ref_col)
        self._alt_col = _human_to_zero(alt_col)

        self._marker_col = _human_to_zero(marker_col)
        self._rsid_col = _human_to_zero(rsid_col)

        self._pval_col = _human_to_zero(pval_col)

        self._beta_col = _human_to_zero(beta_col)
        self._stderr_col = _human_to_zero(stderr_col)

        self._allele_freq_col = _human_to_zero(allele_freq_col)
        self._allele_count_col = _human_to_zero(allele_count_col)
        self._n_samples = _human_to_zero(n_samples_col)

        self._is_log_pval = is_log_pval
        self._is_alt_effect = is_alt_effect

        self.validate_config()

    def validate_config(self):
        """Ensures that a minimally working parser has been created"""
        has_position = (self._marker_col is not None) ^ all(getattr(self, x) is not None
                                                            for x in ('_chr_col', '_pos_col', '_ref_col', '_alt_col'))
        is_valid = has_position and (self._pval_col is not None)
        if not is_valid:
            raise exceptions.ConfigurationException('GWAS parser must specify how to find all required fields')

        if self._allele_count_col and self._allele_freq_col:
            raise exceptions.ConfigurationException('Allele count and frequency options are mutually exclusive')

        if self._allele_count_col and not self._n_samples:
            raise exceptions.ConfigurationException(
                'To calculate allele frequency from counts, you must also provide n_samples')

        return is_valid

    def _split_fields(self, row: str):
        fields = super(GenericGwasLineParser, self)._split_fields(row)
        if len(fields) == 1:
            raise exceptions.LineParseException(
                'Unable to split line into separate fields. This line may have a missing or incorrect delimiter.')
        return fields

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

        # Some optional fields
        beta = None
        stderr_beta = None
        alt_allele_freq = None
        allele_count = None
        n_samples = None

        if self._beta_col is not None:
            beta = values[self._beta_col]

        if self._stderr_col is not None:
            stderr_beta = values[self._stderr_col]

        if self._allele_freq_col is not None:
            alt_allele_freq = values[self._allele_freq_col]

        if self._allele_count_col is not None:
            allele_count = values[self._allele_count_col]
            n_samples = values[self._n_samples]

        # Perform type coercion
        try:
            log_pval = parser_utils.parse_pval_to_log(pval, is_log=self._is_log_pval)
            pos = int(pos)
            if beta is not None:
                beta = None if beta in MISSING_VALUES else float(beta)
            if stderr_beta is not None:
                stderr_beta = None if stderr_beta in MISSING_VALUES else float(stderr_beta)

            if self._allele_freq_col or self._allele_count_col:
                alt_allele_freq = parser_utils.parse_allele_frequency(
                    freq=alt_allele_freq,
                    allele_count=allele_count,
                    n_samples=n_samples,
                    is_alt_effect=self._is_alt_effect
                )
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=values)

        # Some old GWAS files simply won't provide ref or alt information, and the parser will need to do without
        if ref in MISSING_VALUES:
            ref = None

        if alt in MISSING_VALUES:
            alt = None

        return chrom, pos, ref, alt, log_pval, beta, stderr_beta, alt_allele_freq

    def _output_container(self, values):
        return self._container(*values)


class QuickGwasLineParser:
    """
    A "fast path" parser used for pre-standardized input. This parser gains speed by dispensing with
    all the usual re-use and "bad data" tolerance of a more generic tool
    """
    def __init__(self, *, container: ty.Type[BasicVariant] = BasicVariant):
        # The only thing that can be configured is the container (in case we want to support extra fields in the future)
        self._container = container

    @property
    def fields(self) -> ty.Iterable:
        return self._container._fields  # type: ignore

    def __call__(self, row: str) -> BasicVariant:
        # Assume the file format is *exactly* standardized with no extra fields of any kind, no leading or trailing
        #   spaces, and all uses of the delimiter mean what we think they do
        try:
            cols = row.split()
            chrom, pos, ref, alt, log_pvalue = cols[:5]  # These fields are always required

            # For fwd compatibility, the quick-parser will assume that new columns become mandatory & are append-only
            if len(cols) > 5:
                beta = cols[5]
                stderr_beta = cols[6]
                beta = None if beta in MISSING_VALUES else float(beta)
                stderr_beta = None if stderr_beta in MISSING_VALUES else float(stderr_beta)
            else:
                beta = None
                stderr_beta = None

            if len(cols) > 7:
                alt_allele_freq = cols[7]
            else:
                alt_allele_freq = None

            pos = int(pos)
            if ref in MISSING_VALUES:
                ref = None

            if alt in MISSING_VALUES:
                alt = None

            log_pvalue = parser_utils.parse_pval_to_log(log_pvalue, is_log=True)

            alt_allele_freq = parser_utils.parse_allele_frequency(freq=alt_allele_freq, is_alt_effect=True)
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=row)

        return self._container(chrom, pos, ref, alt, log_pvalue, beta, stderr_beta, alt_allele_freq)


####
# Example parsers pre-configured for the LocusZoom standard file format
# Only check the "mandatory" fields
standard_gwas_parser_basic = GenericGwasLineParser(chr_col=1, pos_col=2, ref_col=3, alt_col=4,
                                                   pval_col=5, is_log_pval=True,
                                                   delimiter='\t')

# Parse the "full" standard format (including any additional fields added in the future)
standard_gwas_parser = GenericGwasLineParser(chr_col=1, pos_col=2, ref_col=3, alt_col=4,
                                             pval_col=5, is_log_pval=True,
                                             beta_col=6, stderr_col=7,
                                             allele_freq_col=8,
                                             is_alt_effect=True,
                                             delimiter='\t')

# A "fast" standard parser
standard_gwas_parser_quick = QuickGwasLineParser()
