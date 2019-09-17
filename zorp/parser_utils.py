"""
Utility functions for common parsing or validation operations
"""

import math
import re
import typing as ty

from .const import MISSING_VALUES
from . import exceptions

try:
    from fastnumbers import float
except ImportError:  # pragma: no cover
    pass


REGEX_MARKER = re.compile(r'^(?:chr)?([a-zA-Z0-9]+?):(\d+)[_:]?(\w+)?[/:|]?([^_]+)?_?(.*)?')
REGEX_PVAL = re.compile(r'([\d.\-]+)([\sxeE]*)([0-9\-]*)')


def parse_pval_to_log(value: str, is_neg_log: bool = False) -> ty.Union[float, None]:
    """
    Parse a given number, and return the -log10 pvalue
    """
    if value in MISSING_VALUES or value is None:
        return None

    val = float(value)

    if is_neg_log:  # Take as is
        return val

    # Regular pvalue: validate and convert
    if val < 0 or val > 1:
        raise ValueError('p value is not in the allowed range')

    # 0-values are explicitly allowed and will convert to infinity by design, as they often indicate underflow errors
    #   in the input data.
    if val == 0:
        # Determine whether underflow is due to the source data, or value conversiuon
        if value == '0':
            # The source data is bad, so insert an obvious placeholder value
            return math.inf
        else:
            # h/t @welchr: aggressively turn the underflowing string value into -log10 via regex
            # Only do this if absolutely necessary, because it is a performance hit
            base, _, exponent = REGEX_PVAL.search(value).groups()
            base = float(base)

            if exponent != '':
                exp = float(exponent)
            else:
                exp = 0

            if base == 0:
                return math.inf

            return -(math.log10(float(base)) + float(exp))
    else:
        return -math.log10(val)


def parse_marker(value: str, test: bool = False) -> ty.Union[ty.Tuple[str, str, str, str], None]:
    match = REGEX_MARKER.fullmatch(value)
    if match is not None:
        chrom, pos, ref, alt, _ = match.groups()
        return chrom, pos, ref, alt

    if not test:
        raise exceptions.LineParseException(
            'Could not understand marker format. Must be of format chr:pos or chr:pos_ref/alt')
    else:
        return None


def parse_allele_frequency(*,
                           freq: str = None,
                           allele_count: str = None,
                           n_samples: str = None,
                           is_alt_effect: bool = True) -> ty.Union[float, None]:
    """
    Parse an allele frequency, OR convert counts to frequency.
    :param freq:
    :param allele_count:
    :param n_samples:
    :param is_alt_effect:
    :return:
    """
    if freq is not None and allele_count is not None:
        # A single uber-func is generally less than the performance penalty of calling two separate functions
        raise exceptions.ConfigurationException('Frequency and allele count options are mutually exclusive')

    if freq is None and (allele_count in MISSING_VALUES or n_samples in MISSING_VALUES):  # Allele count parsing
        return None
    elif freq is None and allele_count is not None:
        allele_freq = int(allele_count) / int(n_samples) / 2  # 2 alleles per sample
    elif freq in MISSING_VALUES:  # Frequency-based parsing
        return None
    else:
        allele_freq = float(freq)

    # No matter how the frequency is specified, this stuff is always done
    if allele_freq < 0 or allele_freq > 1:
        raise ValueError('Allele frequency is not in the allowed range')

    if not is_alt_effect:  # Orient the frequency to the alt allele
        return 1 - allele_freq
    else:
        return allele_freq


def human_to_zero(value):
    """
    Many of our parsers use 1-based indices as arguments, and convert to 0-based indices for internal usage
        (eg, "column 1 = index 0 in this list of fields")
    """
    if value is None:
        return value
    else:
        return value - 1
