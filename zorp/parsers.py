"""
Parsers: handle the act of reading one entity (such as line)
"""
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
        self.chrom = chrom  # type: str
        self.pos = pos  # type: int
        self.rsid = rsid  # type: str
        self.ref = ref  # type: str
        self.alt = alt  # type: str
        self.neg_log_pvalue = neg_log_pvalue  # type: float

        # # Optional fields for future expansion
        # af: float
        self.beta = beta  # type: float
        self.stderr_beta = stderr_beta  # type: float

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
        ref_alt = '_{}/{}'.format(self.ref, self.alt) \
            if (self.ref and self.alt) else ''
        return '{}:{}{}'.format(self.chrom, self.pos, ref_alt)

    def to_dict(self):
        # Some tools expect the data in a mutable form (eg dicts)
        return {s: getattr(self, s, None) for s in self._fields}


def TupleLineParser(*args, container: ty.Callable = tuple, delimiter='\t', **kwargs):
    """
    Parse a line of text and return a tuple of the fields. Performs no type coercion

    This isn't recommended for everyday parsing, but it is useful internally for format detection (where we need to
        split columns of data, but aren't yet ready to clean and assign meaning to the values)
    """
    def inner(line: str):
        """Return a stateful closure that actually does the work of parsing"""
        try:
            values = line.strip().split(delimiter)
            return container(values)
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=line)
    return inner


def GenericGwasLineParser(
        *args,
        delimiter: str = '\t',
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
    """
    A simple parser that extracts GWAS information from a flexible file format.

    Constructor expects human-friendly column numbers (first = column 1)
    """
    def validate_config():
        """Ensures that a minimally working parser has been created"""
        # Some old gwas files may not have ref and alt (incomplete marker, or missing columns). These fields aren't
        #   strictly required, but we really really like to have them
        has_position = (_marker_col is not None) ^ all(x is not None
                                                       for x in (_chrom_col, _pos_col))
        # If we do have one allele, we must have both
        both_markers = (_ref_col is None and _alt_col is None) or \
                       (_ref_col is not None and _alt_col is not None)

        is_valid = has_position and both_markers and (_pvalue_col is not None)
        if not is_valid:
            raise exceptions.ConfigurationException('GWAS parser must specify how to find all required fields')

        if _allele_count_col is not None and _allele_freq_col is not None:
            raise exceptions.ConfigurationException('Allele count and frequency options are mutually exclusive')

        if _allele_count_col is not None and _n_samples_col is None:
            raise exceptions.ConfigurationException(
                'To calculate allele frequency from counts, you must also provide n_samples')

        return is_valid

    def inner(line):
        # Return a stateful closure that does the actual work of parsing
        try:
            fields = line.strip().split(delimiter)
            if len(fields) == 1:
                raise exceptions.LineParseException(
                    'Unable to split line into separate fields. This line may have a missing or incorrect delimiter.')

            # Fetch values
            ref = None
            alt = None
            if _marker_col is not None:
                chrom, pos, ref, alt = utils.parse_marker(fields[_marker_col])
            else:
                chrom = fields[_chrom_col]
                pos = fields[_pos_col]

            if chrom.startswith('chr'):
                chrom = chrom[3:]

            chrom = chrom.upper()

            # Explicit columns will override a value from the marker, by design
            if _ref_col is not None:
                ref = fields[_ref_col]

            if _alt_col is not None:
                alt = fields[_alt_col]

            pval = fields[_pvalue_col]

            # Some optional fields
            rsid = None
            beta = None
            stderr_beta = None
            alt_allele_freq = None
            allele_count = None
            n_samples = None

            if _rsid_col is not None:
                rsid = fields[_rsid_col]
                if rsid in MISSING_VALUES:
                    rsid = None
                elif not rsid.startswith('rs'):
                    rsid = 'rs' + rsid

            if _beta_col is not None:
                beta = fields[_beta_col]

            if _stderr_col is not None:
                stderr_beta = fields[_stderr_col]

            if _allele_freq_col is not None:
                alt_allele_freq = fields[_allele_freq_col]

            if _allele_count_col is not None:
                allele_count = fields[_allele_count_col]
                n_samples = fields[_n_samples_col]

            # Perform type coercion
            log_pval = utils.parse_pval_to_log(pval, is_neg_log=_is_neg_log_pvalue)

            try:
                pos = int(pos)
            except ValueError:
                # Some programs seem to write long positions using scientific notation, which int cannot handle
                try:
                    pos = int(float(pos))
                except ValueError:
                    # If we still can't parse, it's probably bad data
                    raise exceptions.LineParseException(
                        'Positions should be specified as integers. Could not parse value: {}'.format(pos))

            if beta is not None:
                beta = None if beta in MISSING_VALUES else float(beta)
            if stderr_beta is not None:
                stderr_beta = None if stderr_beta in MISSING_VALUES else float(stderr_beta)

            if _allele_freq_col or _allele_count_col:
                alt_allele_freq = utils.parse_allele_frequency(
                    freq=alt_allele_freq,
                    allele_count=allele_count,
                    n_samples=n_samples,
                    is_alt_effect=_is_alt_effect
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

            result = container(chrom, pos, rsid, ref, alt, log_pval, beta, stderr_beta, alt_allele_freq)
        except Exception as e:
            raise exceptions.LineParseException(str(e), line=line)
        return result

    # Convert the user-provided values to field array indices (0-based), and validate config
    # Chrom field has a legacy alias, allowing older parser configs to work.
    inner._chrom_col = _chrom_col = utils.human_to_zero(chrom_col) if chrom_col is not None else utils.human_to_zero(chr_col)  # type: ignore  # noqa
    inner._pos_col = _pos_col = utils.human_to_zero(pos_col)  # type: ignore
    inner._ref_col = _ref_col = utils.human_to_zero(ref_col)  # type: ignore
    inner._alt_col = _alt_col = utils.human_to_zero(alt_col)  # type: ignore

    inner._marker_col = _marker_col = utils.human_to_zero(marker_col)  # type: ignore

    # Support legacy alias for field name
    inner._pvalue_col = _pvalue_col = utils.human_to_zero(pvalue_col) if pvalue_col is not None else utils.human_to_zero(pval_col)  # type: ignore  # noqa

    inner._rsid_col = _rsid_col = utils.human_to_zero(rsid_col)  # type: ignore
    inner._beta_col = _beta_col = utils.human_to_zero(beta_col)  # type: ignore
    inner._stderr_col = _stderr_col = utils.human_to_zero(stderr_beta_col)  # type: ignore

    inner._allele_freq_col = _allele_freq_col = utils.human_to_zero(allele_freq_col)  # type: ignore
    inner._allele_count_col = _allele_count_col = utils.human_to_zero(allele_count_col)  # type: ignore
    inner._n_samples_col = _n_samples_col = utils.human_to_zero(n_samples_col)  # type: ignore

    # The latter option is an alias for legacy reasons
    inner._is_neg_log_pvalue = _is_neg_log_pvalue = is_neg_log_pvalue or is_log_pval  # type: ignore
    inner._is_alt_effect = _is_alt_effect = is_alt_effect  # type: ignore

    # Raise an exception if the provided options are invalid
    validate_config()

    # Provide the outside world with access to additional named attributes
    # We are slightly abusing closures, but the end result is ~10% is faster than a class-based callable
    inner.fields = container._fields  # type: ignore
    return inner
