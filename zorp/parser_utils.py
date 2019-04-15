"""
Utility functions for common parsing or validation operations
"""

import math
import re

from . import exceptions


REGEX_MARKER = re.compile(r'(?:chr)?(.+):(\d+)[_:]?(\w+)?[/:|]?([^_]+)?_?(.*)?')


def parse_pval_to_log(value, is_log=False):
    """Parse a given number, and return the -log10 pvalue"""
    if value is None:
        return value
    val = float(value)  # TODO: use Ryan's utility to handle exponents

    if is_log:  # Take as is
        return val

    # Regular pvalue: validate and convert
    if val < 0 or val > 1:
        raise ValueError('p value is not in the allowed range')

    # 0-values are explicitly allowed and will convert to infinity by design
    if val == 0:
        return math.inf
    else:
        return -math.log10(val)


def parse_marker(value: str, test: bool = False):
    match = REGEX_MARKER.fullmatch(value)
    if match:
        return match.groups()

    if not test:
        raise exceptions.LineParseException('Could not understand marker format. Must be of format chr:pos or chr:pos_ref/alt')
    else:
        return None


