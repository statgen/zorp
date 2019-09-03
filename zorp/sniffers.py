"""
Auto-detect the appropriate options for parsing a GWAS file
"""

# pvalue or logpvalue
# How many header rows? (in a text file, we don't necessarily know this)
# What are the headers names?
# Get columns for chrom/pos/ref/alt
import binascii
import itertools
import typing as ty

try:
    from fastnumbers import float
except ImportError:
    pass

from .const import MISSING_VALUES
from . import (
    const,
    exceptions,
    parsers,
    parser_utils,
    readers
)


def is_numeric(val: str) -> bool:
    """Check whether an unparsed string is a numeric value"""
    if val in MISSING_VALUES:
        return True

    try:
        float(val)
    except Exception:
        return False
    else:
        return True


def is_header(row: str, *, comment_char="#", delimiter='\t') -> bool:
    """
    This assumes two basic rules: the line is not a comment, and gwas data is more likely to be numeric than headers
    """
    return row.startswith(comment_char) or all(not is_numeric(field)
                                               for field in row.split(delimiter))


def levenshtein(s1, s2):
    # CC_BY_SA https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = list(range(len(s2) + 1))
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            # j+1 instead of j since previous_row and current_row are one character longer than s2
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]


def find_column(column_synonyms: tuple, header_names: list, threshold: int = 2):
    #  Find the column name that best matches
    best_score = threshold + 1
    best_match = None
    for i, header in enumerate(header_names):
        if header is None:
            # If header is empty, don't consider it for a match
            # Nulling a header provides a way to exclude something from future searching
            continue

        score = min(levenshtein(header, s) for s in column_synonyms)
        if score < best_score:
            best_score = score
            best_match = i
    return best_match


def get_pval_column(header_names: list, data_rows: ty.Iterable) \
        -> ty.Union[dict, None]:
    """
    Return (column, is_log), or None/None if config not found
    Column numbers are 1-based, for compatibility with parsers
    """
    LOGPVALUE_FIELDS = ('neg_log_pvalue', 'log_pvalue', 'log_pval', 'logpvalue')
    PVALUE_FIELDS = ('pvalue', 'p.value', 'pval', 'p_score', 'p')

    data = itertools.islice(data_rows, 100)

    def _validate_p(col: int, data: ty.Iterator, is_log: bool) -> bool:
        # All values must be parseable
        vals = [row[col] for row in data]
        cleaned_vals = [None if val in const.MISSING_VALUES else val
                        for val in vals]
        try:
            for v in cleaned_vals:
                parser_utils.parse_pval_to_log(v, is_log=is_log)
        except Exception:
            return False

        return True

    log_p_col = find_column(LOGPVALUE_FIELDS, header_names)
    p_col = find_column(PVALUE_FIELDS, header_names)

    if log_p_col is not None and _validate_p(log_p_col, data, True):
        return {'pval_col': log_p_col + 1, 'is_log_pval': True}
    elif p_col is not None and _validate_p(p_col, data, False):
        return {'pval_col': p_col + 1, 'is_log_pval': False}

    # Could not auto-determine an appropriate pvalue column
    return None


def get_chrom_pos_ref_alt_columns(header_names: list, data_rows: ty.Iterable):
    """
    Find the information required to uniquely identify the marker. Returns 1-based column indices for compatibility
        with line parsers.
    :param header_names:
    :param data_rows:
    :return:
    """
    # Get from either a marker, or 4 separate columns
    MARKER_FIELDS = ('snpid', 'marker', 'markerid', 'snpmarker', 'chr:position')
    CHR_FIELDS = ('chrom', 'chr')
    POS_FIELDS = ('position', 'pos', 'begin', 'beg', 'bp', 'end', 'ps')

    data = itertools.islice(data_rows, 100)

    # TODO: How to handle orienting ref vs effect?
    # Order matters: consider ambiguous field names for ref before alt
    REF_FIELDS = ('A1', 'ref', 'reference', 'allele0', 'allele1')
    ALT_FIELDS = ('A2', 'alt', 'alternate', 'allele1', 'allele2')

    first_row = next(data)
    marker_col = find_column(MARKER_FIELDS, header_names)

    if marker_col is not None and parser_utils.parse_marker(first_row[marker_col], test=True):
        return {'marker_col': marker_col + 1}

    # If single columns were incomplete, attempt to auto detect 4 separate columns. All 4 must
    #  be found for this function to report a match.
    headers_marked = header_names.copy()
    to_find = [
        ['chr_col', CHR_FIELDS],
        ['pos_col', POS_FIELDS],
        ['ref_col', REF_FIELDS],
        ['alt_col', ALT_FIELDS],
    ]
    config = {}
    for col_name, col_choices in to_find:
        col = find_column(col_choices, headers_marked)
        if col is None:
            return None

        config[col_name] = col + 1
        # Once a column has been assigned, remove it from consideration for future matches
        headers_marked[col] = None

    return config


def get_effect_size_columns(header_names: list, data_rows: ty.Iterable):
    BETA_FIELDS = ('beta', 'effect_size', 'alt_effsize', 'effect')
    STDERR_BETA_FIELDS = ('stderr_beta', 'stderr', 'effect_size_sd', 'se')

    data = itertools.islice(data_rows, 100)

    def _validate_numeric(col: int, data: ty.Iterator) -> bool:
        vals = [row[col] for row in data]
        cleaned_vals = [val for val in vals
                        if val not in const.MISSING_VALUES]
        try:

            for v in cleaned_vals:
                float(v)
        except:
            return False
        return True

    beta_col = find_column(BETA_FIELDS, header_names, threshold=0)
    stderr_col = find_column(STDERR_BETA_FIELDS, header_names, threshold=0)

    ret = {}
    if beta_col is not None and _validate_numeric(beta_col, data):
        ret['beta_col'] = beta_col + 1

    if stderr_col is not None and _validate_numeric(stderr_col, data):
        ret['stderr_col'] = stderr_col + 1

    return ret or None


def get_reader(filename: ty.Union[ty.Iterable, str]) -> ty.Type[readers.BaseReader]:
    """Suggest an appropriate reader class based on source of data: iterable, tabix, or text file"""
    if not isinstance(filename, str):
        return readers.IterableReader

    with open(filename, 'rb') as test_f:
        is_gz = binascii.hexlify(test_f.read(2)) == b'1f8b'

    if is_gz:
        return readers.TabixReader
    else:
        return readers.TextFileReader


def get_headers(reader, comment_char: str = "#", delimiter: str = '\t', max_check=100) \
        -> ty.Tuple[int, ty.Union[str, None]]:
    """Identify the number of header rows, and the content of the one likely to contain column headers"""
    last_row = None
    for i, row in enumerate(reader):
        # TODO: Move some of this method to parser class to keep the reader domain-agnostic
        if not is_header(row, comment_char=comment_char, delimiter=delimiter):
            return i, last_row
        elif i > max_check:
            raise exceptions.SnifferException(f'No headers found after limit of {max_check} rows')
        last_row = row

    raise exceptions.SnifferException('No headers found after searching entire file')


def guess_gwas(filename: ty.Union[ty.Iterable, str], *,
               skip_rows=None,
               parser: parsers.AbstractLineParser = None,
               delimiter: str = '\t',
               **kwargs) -> readers.BaseReader:
    """
    Read tab delimited rows of data and return a fully configured reader/parser
    This will try to identify the following options:
    - Type of file (text, gzip, etc)
    - Number of header rows (unless provided explicitly)
    - How to parse the file (unless a parser is explicitly provided)

    Supports receiving an iterable (instead of filename), primarily to support unit testing
    """
    reader_class = get_reader(filename)
    n_headers, header_text = get_headers(reader_class(filename, parser=None), delimiter=delimiter)

    # Don't try to guess parser options if options are explicitly provided. That would be silly.
    to_skip = n_headers if skip_rows is None else skip_rows
    if parser is None:
        # FIXME: Handle case of files with no header rows
        header_names = header_text.lower().lstrip('#').split(delimiter)

        # The first effort at field detection just extracts fields, with no value cleanup
        parser = parsers.TupleLineParser(delimiter=delimiter)

        # Any kwargs not specified for this function are assumed to be reader options, and passed directly in
        data_reader = reader_class(filename, skip_rows=to_skip, parser=parser, **kwargs)

        p_config = get_pval_column(header_names, data_reader)
        if p_config is None:
            raise exceptions.SnifferException('Could not find required field: pvalue')

        header_names[p_config['pval_col'] - 1] = None  # Remove this column from consideration for other matches
        position_config = get_chrom_pos_ref_alt_columns(header_names, data_reader)

        if position_config is None:
            raise exceptions.SnifferException('Could not find SNP identifier columns (position or marker)')

        for v in position_config.values():
            header_names[v - 1] = None  # Remove columns from consideration

        beta_config = get_effect_size_columns(header_names, data_reader)

        # Configure a reader and parser based on the auto-detected file options
        options = {**p_config, **position_config, **(beta_config or {})}

        parser = parsers.GenericGwasLineParser(**options)

    return reader_class(filename, skip_rows=to_skip, parser=parser)
