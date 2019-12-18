"""
Parsers: handle the act of reading one entity (such as line)
"""
import abc
import math
import numbers
import typing as ty

try:
    from fastnumbers import int, float
except ImportError:  # pragma: no cover
    pass

from .const import MISSING_VALUES
from . import exceptions, parser_utils as utils


class BasicVariant:
    """
    Store GWAS results in a predictable format, with a minimal set of fields; optimize for name-based attribute access
    """
    # Slots specify the data  this holds (a performance optimization); _fields is human-curated list
    __slots__ = ('chrom', 'pos', 'ref', 'alt', 'neg_log_pvalue', 'beta', 'stderr_beta', 'alt_allele_freq', 'rsid')
    _fields = ('chrom', 'pos', 'rsid', 'ref', 'alt', 'neg_log_pvalue', 'beta', 'stderr_beta', 'alt_allele_freq')

    def __init__(self, chrom, pos, rsid, ref, alt, neg_log_pvalue, beta, stderr_beta, alt_allele_freq):
        self.chrom: str = chrom
        self.pos: int = pos
        self.rsid: str = rsid
        self.ref: str = ref
        self.alt: str = alt
        self.neg_log_pvalue: float = neg_log_pvalue

        # # Optional fields for future expansion
        # af: float
        self.beta: float = beta
        self.stderr_beta: float = stderr_beta

        self.alt_allele_freq = alt_allele_freq

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
    def maf(self) -> ty.Union[numbers.Number, None]:
        af = self.alt_allele_freq
        return min(af, 1 - af) if af is not None else None

    @property
    def marker(self) -> str:
        """Specify the marker in a string format compatible with UM LD server and other variant-specific requests"""
        ref_alt = f'_{self.ref}/{self.alt}' if (self.ref and self.alt) else ''
        return f'{self.chrom}:{self.pos}{ref_alt}'

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
        :param container: A data structure (eg tuple or slotted object) that will be populated with the parsed results
        """
        self._container = container

    @abc.abstractmethod
    def __call__(self, row: str) -> object:
        pass


class TupleLineParser(AbstractLineParser):
    """
    Parse a line of text and return a tuple of the fields. Performs no type coercion

    This isn't recommended for everyday parsing, but it is useful internally for format detection (where we need to
        split columns of data, but aren't yet ready to clean and assign meaning to the values)
    """
    def __init__(self, *args, container: ty.Callable = tuple, delimiter='\t', **kwargs):
        super(TupleLineParser, self).__init__(*args, container=container, **kwargs)
        self._delimiter = delimiter

    def __call__(self, line: str):
        try:
            values = line.strip().split(self._delimiter)
            return self._container(values)
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=line)


class GenericGwasLineParser(TupleLineParser):
    """
    A simple parser that extracts GWAS information from a flexible file format.

    Constructor expects human-friendly column numbers (first = column 1)
    """
    # Provide faster attribute access because these field lookups will happen a *lot*
    __slots__ = [
        '_container',
        '_chrom_col', 'rsid_col', '_pos_col', '_ref_col', '_alt_col', '_marker_col', '_pvalue_col',
        '_pvalue_col', '_beta_col', '_stderr_col', '_allele_freq_col', '_allele_count_col',
        '_n_samples_col', '_is_neg_log_pvalue', '_is_alt_effect'
    ]

    def __init__(self, *args,
                 container: ty.Callable[..., BasicVariant] = BasicVariant,
                 # Variant identifiers: marker OR individual
                 chrom_col: int = None, chr_col: int = None,  # Legacy alias
                 pos_col: int = None,
                 ref_col: int = None,
                 alt_col: int = None,
                 marker_col: int = None,
                 # Other required data
                 pvalue_col: int = None, pval_col: int = None,  # Legacy alias
                 # Optional fields
                 rsid_col: int = None,
                 beta_col: int = None,
                 stderr_beta_col: int = None,

                 allele_freq_col: int = None,  # As freq OR by two count fields
                 allele_count_col: int = None,
                 n_samples_col: int = None,
                 # Other configuration options that apply to every row as constants
                 is_neg_log_pvalue: bool = False, is_log_pval: bool = False,  # Legacy alias
                 is_alt_effect: bool = True,  # whether effect allele is oriented towards alt
                 **kwargs):
        super().__init__(*args, **kwargs)

        # All GWAS parsers can specify
        self._container = container

        # Chrom field has a legacy alias, allowing older parser configs to work.
        self._chrom_col = utils.human_to_zero(chrom_col) if chrom_col is not None else utils.human_to_zero(chr_col)
        self._pos_col = utils.human_to_zero(pos_col)
        self._ref_col = utils.human_to_zero(ref_col)
        self._alt_col = utils.human_to_zero(alt_col)

        self._marker_col = utils.human_to_zero(marker_col)

        # Support legacy alias for field name
        self._pvalue_col = utils.human_to_zero(pvalue_col) if pvalue_col is not None else utils.human_to_zero(pval_col)

        self._rsid_col = utils.human_to_zero(rsid_col)
        self._beta_col = utils.human_to_zero(beta_col)
        self._stderr_col = utils.human_to_zero(stderr_beta_col)

        self._allele_freq_col = utils.human_to_zero(allele_freq_col)
        self._allele_count_col = utils.human_to_zero(allele_count_col)
        self._n_samples_col = utils.human_to_zero(n_samples_col)

        self._is_neg_log_pvalue = is_neg_log_pvalue or is_log_pval  # The latter option is an alias for legacy reasons
        self._is_alt_effect = is_alt_effect

        self.validate_config()

    def validate_config(self):
        """Ensures that a minimally working parser has been created"""
        # Some old gwas files may not have ref and alt (incomplete marker, or missing columns). These fields aren't
        #   strictly required, but we really really like to have them
        has_position = (self._marker_col is not None) ^ all(getattr(self, x) is not None
                                                            for x in ('_chrom_col', '_pos_col'))
        # If we do have one allele, we must have both
        both_markers = (self._ref_col is None and self._alt_col is None) or \
                       (self._ref_col is not None and self._alt_col is not None)

        is_valid = has_position and both_markers and (self._pvalue_col is not None)
        if not is_valid:
            raise exceptions.ConfigurationException('GWAS parser must specify how to find all required fields')

        if self._allele_count_col and self._allele_freq_col:
            raise exceptions.ConfigurationException('Allele count and frequency options are mutually exclusive')

        if self._allele_count_col and not self._n_samples_col:
            raise exceptions.ConfigurationException(
                'To calculate allele frequency from counts, you must also provide n_samples')

        return is_valid

    @property
    def fields(self) -> ty.Container:
        return self._container._fields  # type: ignore

    def __call__(self, line: str):
        try:
            fields = line.strip().split(self._delimiter)
            if len(fields) == 1:
                raise exceptions.LineParseException(
                    'Unable to split line into separate fields. This line may have a missing or incorrect delimiter.')

            # Fetch values
            ref = None
            alt = None
            if self._marker_col is not None:
                chrom, pos, ref, alt = utils.parse_marker(fields[self._marker_col])
            else:
                chrom = fields[self._chrom_col]
                pos = fields[self._pos_col]

            if chrom.startswith('chr'):
                chrom = chrom[3:]

            chrom = chrom.upper()

            # Explicit columns will override a value from the marker, by design
            if self._ref_col is not None:
                ref = fields[self._ref_col]

            if self._alt_col is not None:
                alt = fields[self._alt_col]

            pval = fields[self._pvalue_col]

            # Some optional fields
            rsid = None
            beta = None
            stderr_beta = None
            alt_allele_freq = None
            allele_count = None
            n_samples = None

            if self._rsid_col is not None:
                rsid = fields[self._rsid_col]
                if rsid in MISSING_VALUES:
                    rsid = None
                elif not rsid.startswith('rs'):
                    rsid = 'rs' + rsid

            if self._beta_col is not None:
                beta = fields[self._beta_col]

            if self._stderr_col is not None:
                stderr_beta = fields[self._stderr_col]

            if self._allele_freq_col is not None:
                alt_allele_freq = fields[self._allele_freq_col]

            if self._allele_count_col is not None:
                allele_count = fields[self._allele_count_col]
                n_samples = fields[self._n_samples_col]

            # Perform type coercion
            log_pval = utils.parse_pval_to_log(pval, is_neg_log=self._is_neg_log_pvalue)

            try:
                pos = int(pos)
            except ValueError:
                # Some programs seem to write long positions using scientific notation, which int cannot handle
                try:
                    pos = int(float(pos))
                except ValueError:
                    # If we still can't parse, it's probably bad data
                    raise exceptions.LineParseException(
                        f'Positions should be specified as integers. Could not parse value: {pos}')

            if beta is not None:
                beta = None if beta in MISSING_VALUES else float(beta)
            if stderr_beta is not None:
                stderr_beta = None if stderr_beta in MISSING_VALUES else float(stderr_beta)

            if self._allele_freq_col or self._allele_count_col:
                alt_allele_freq = utils.parse_allele_frequency(
                    freq=alt_allele_freq,
                    allele_count=allele_count,
                    n_samples=n_samples,
                    is_alt_effect=self._is_alt_effect
                )

            # Some old GWAS files simply won't provide ref or alt information, and the parser will need to do without
            if ref in MISSING_VALUES:
                ref = None

            if isinstance(ref, str):
                ref = ref.upper()

            if alt in MISSING_VALUES:
                alt = None

            if isinstance(alt, str):
                alt = alt.upper()

            container = self._container(chrom, pos, rsid, ref, alt, log_pval, beta, stderr_beta, alt_allele_freq)
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=line)
        return container
