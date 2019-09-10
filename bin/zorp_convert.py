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
    sniffers,
)

logger = logging.getLogger(__name__)


def main(source: ty.Union[str, ty.Iterable],
         out_fn: str,
         parser_options: dict,
         auto_config=False,
         skip_rows=None,
         skip_errors=True,
         max_errors=100,
         make_tabix: bool = False) -> str:
    try:
        parser = parsers.GenericGwasLineParser(**parser_options)
    except exceptions.ConfigurationException:
        parser = None

    if source is None:
        source = sys.stdin

    if not auto_config and (skip_rows is None or parser is None):
        logger.error('Please provide all options required to parse the file, or use the --auto flag to guess')
        sys.exit(1)

    reader = sniffers.guess_gwas_generic(source, skip_rows=skip_rows, parser=parser,
                                         skip_errors=skip_errors, max_errors=max_errors)

    # By design, this filters all missing pvalues out of the result file
    reader.add_filter('log_pvalue', lambda v, _: v is not None)

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
    parser.add_argument('--skip-rows', dest='skip_rows', type=int, default=None,
                        help='The number of non-data lines that will be skipped')
    parser.add_argument('--max-errors', type=int, default=100, dest='max_errors',
                        help='How many unparseable lines can be skipped before the script stops')
    parser.add_argument('--stop-on-error', action='store_true', dest='stop_on_error',
                        help='Should we stop parsing immediately on the first error?')

    # Offering a filtering option might be useful, but hold this feature until it's more useful than tabix
    # parser.add_argument('--filter', action='append', nargs=2,
    #                     help='Only accept lines that match the provided string value')

    # Format specification options  # TODO: Improve argparse to make clear which options are mutually exclusive.
    # There are multiple ways of specifying chrom/pos/ref/alt columns: auto-detect, marker, or explicit columns
    parser.add_argument('--auto', action='store_true', dest='auto',
                        help='Automatically detect appropriate parser options')

    parser.add_argument('-m', '--marker_col', type=int, dest='marker_col',
                        help='Column number with marker (used to find chr, pos, ref, and alt)')

    parser.add_argument('-c', '--chrom_col', type=int, dest='chrom_col',
                        help='Column number with chromosome')
    parser.add_argument('-b', '--pos_col', type=int, dest='pos_col',
                        help='Column number with position')
    parser.add_argument('-r', '--ref_col', type=int, dest='ref_col',
                        help='Column number with refvar')
    parser.add_argument('-a', '--alt_col', type=int, dest='alt_col',
                        help='Column number with altvar')

    # Optional fields
    parser.add_argument('--beta_col', type=int, dest='beta_col',
                        help='Column number with beta (effect size)')
    parser.add_argument('--stderr_beta_col', type=int, dest='stderr_beta_col',
                        help='Column number with stderr of effect size')

    # This argument is always required
    parser.add_argument('-p', '--pvalue_col', type=int, dest='pvalue_col',
                        help='Column number with pvalue (or logpvalue)')
    parser.add_argument('--is-neg-log-pvalue', action='store_true', dest='is_neg_log_pvalue',
                        help='Is the pvalue data stored in -log10(p) form?')

    parser.add_argument('--compress', action='store_true', dest='compress',
                        help='Whether to tabix and compress the output')

    args = parser.parse_args()

    main(
        args.file,
        args.dest,
        {
            'marker_col': args.marker_col,
            'chrom_col': args.chrom_col,
            'pos_col': args.pos_col,
            'ref_col': args.ref_col,
            'alt_col': args.alt_col,
            'pvalue_col': args.pvalue_col,
            'is_neg_log_pvalue': args.is_neg_log_pvalue,

            # Optional args
            'beta_col': args.beta_col,
            'stderr_beta_col': args.stderr_beta_col
        },
        auto_config=args.auto,
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
