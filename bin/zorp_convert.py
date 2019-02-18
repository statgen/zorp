#!/usr/bin/env python3
"""
Command line utility that reads a GWAS file and converts it to a standardized output

Supports tabix files (via command line argument), or any text stream (eg piping from another program)
"""

import argparse
import logging
import sys
import typing as ty

from zorp import (
    exceptions,
    parsers,
    readers
)

logger = logging.getLogger(__name__)


def main(source: ty.Union[str, ty.Iterable],
         out_fn: str,
         parser_options: dict,
         skip_rows=0,
         skip_errors=True,
         max_errors=100,
         make_tabix: bool = False) -> str:
    parser = parsers.GenericGwasLineParser(**parser_options)

    # TODO: Add support for file format detection- not every infile will be tabix/gz!
    if source is None:
        reader_class = readers.IterableReader
        source = sys.stdin
    else:
        reader_class = readers.TabixReader

    reader = reader_class(source,
                          parser=parser,
                          skip_errors=skip_errors, skip_rows=skip_rows, max_errors=max_errors)

    reader.add_filter('pvalue', lambda v, _: v is not None)

    try:
        dest_fn = reader.write(out_fn, make_tabix=make_tabix)
    except exceptions.TooManyBadLinesException:
        logger.error('ERROR: Too many lines failed to parse; stopping.')
    except Exception:
        logger.exception('Conversion failed due to unknown error')
    else:
        logger.info('Conversion succeeded! Results written to: {}'.format(dest_fn))
        return dest_fn
    finally:
        for n, reason, _ in reader.errors:
            logger.error('Excluded row {} from output due to parse error: {}'.format(n, reason))


def run():
    """Command line arguments"""
    parser = argparse.ArgumentParser(description="A utility to convert GWAS files into a standardized format."
                                                 "All column numbers should specify human-friendly indices (1-based)")
    parser.add_argument('file', nargs='?', default=None,
                        help='The path to the input GWAS file')
    parser.add_argument('--dest', type=str, required=True,
                        help='The filename to use when saving the output GWAS file')
    parser.add_argument('--skip-rows', dest='skip_rows', type=int, default=0,
                        help='The number of non-data lines that will be skipped')
    parser.add_argument('--max-errors', type=int, default=100, dest='max_errors',
                        help='How many unparseable lines can be skipped before the script stops')
    parser.add_argument('--stop-on-error', action='store_true', dest='stop_on_error',
                        help='Should we stop parsing immediately on the first error?')

    # Format specification options
    parser.add_argument('-c', '--chr_col', type=int, dest='chr_col', required=True,
                        help='Column number with chromosome')
    parser.add_argument('-b', '--pos_col', type=int, dest='pos_col', required=True,
                        help='Column number with position')
    parser.add_argument('-r', '--ref_col', type=int, dest='ref_col', required=True,
                        help='Column number with refvar')
    parser.add_argument('-a', '--alt_col', type=int, dest='alt_col', required=True,
                        help='Column number with altvar')
    parser.add_argument('-p', '--pval_col', type=int, dest='pval_col', required=True,
                        help='Column number with pvalue (or logpvalue)')
    parser.add_argument('--is-log-pval', action='store_true', dest='is_log_pval',
                        help='Is the pvalue data stored in -log10(p) form?')

    parser.add_argument('--compress', action='store_true', dest='compress',
                        help='Whether to tabix and compress the output')

    args = parser.parse_args()

    main(
        args.file,
        args.dest,
        {
            'chr_col': args.chr_col,
            'pos_col': args.pos_col,
            'ref_col': args.ref_col,
            'alt_col': args.alt_col,
            'pval_col': args.pval_col,
            'is_log_pval': args.is_log_pval,
        },
        max_errors=args.max_errors,
        skip_errors=not args.stop_on_error,
        skip_rows=args.skip_rows,
        make_tabix=args.compress
    )


if __name__ == '__main__':
    import time

    t1 = time.time()
    run()
    logger.info('Conversion complete: Elapsed time: {}s'.format(time.time() - t1))
