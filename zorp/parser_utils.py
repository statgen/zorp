"""
Utility functions for common parsing or validation operations
"""

import math
import re

from .const import MISSING_VALUES
from . import exceptions

try:
    from fastnumbers import float
except ImportError:
    pass


REGEX_MARKER = re.compile(r'^(?:chr)?([a-zA-Z0-9]+?):(\d+)[_:]?(\w+)?[/:|]?([^_]+)?_?(.*)?')
REGEX_PVAL = re.compile("([\d.\-]+)([\sxeE]*)([0-9\-]*)")


def parse_pval_to_log(value, is_log=False):
    """
    Parse a given number, and return the -log10 pvalue
    `is_log` should really be "is negative log", and is confusingly named for legacy reasons. FIXME: Change that
    """
    if value in MISSING_VALUES or value is None:
        return None

    val = float(value)

    if is_log:  # Take as is
        return val

    # Regular pvalue: validate and convert
    if val < 0 or val > 1:
        raise ValueError('p value is not in the allowed range')

    # 0-values are explicitly allowed and will convert to infinity by design, as they often indicate underflow errors
    #   in the input data.
    if val == 0:
        # Determine whether underflow is due to the source data, or due to python reading in the number
        if value == '0':
            # The source data is bad, so insert an obvious placeholder value
            return math.inf
        else:
            # h/t @welchr: aggressively turn the underflowing string value into -log10 via regex
            # Only do this if absolutely necessary, because it is a performance hit
            base, _, exponent = REGEX_PVAL.search(value).groups()
            base = float(base)

            if exponent != '':
                exponent = float(exponent)
            else:
                exponent = 0

            if base == 0:
                return math.inf

            return -(math.log10(float(base)) + float(exponent))
    else:
        return -math.log10(val)


def parse_marker(value: str, test: bool = False):
    match = REGEX_MARKER.fullmatch(value)
    if match:
        chrom, pos, ref, alt, _ = match.groups()
        return chrom, pos, ref, alt

    if not test:
        raise exceptions.LineParseException(
            'Could not understand marker format. Must be of format chr:pos or chr:pos_ref/alt')
    else:
        return None
