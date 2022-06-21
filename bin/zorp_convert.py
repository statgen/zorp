#!/usr/bin/env python3
"""
Command line utility that reads a GWAS file and converts it to a standardized output

Supports tabix files (via command line argument), or any text stream (eg piping from another program)
"""

import argparse
import functools
import logging
import sys
import traceback
import typing as ty

from zorp import (
    exceptions,
    parsers,
    sniffers,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(message)s')


def suppress_broken_pipe_msg(f):
    """
    When using the CLI as part of a script pipeline, we want to gracefully handle the next command closing the pipe
     early (eg head). This is a workaround for the fact that python prints an error message to the console even
     when an error is correctly handled.
    https://stackoverflow.com/questions/14207708/ioerror-errno-32-broken-pipe-python
    """
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except SystemExit:
            raise
        except Exception:
            traceback.print_exc()
            exit(1)
        finally:
            try:
                sys.stdout.flush()
            finally:
                try:
                    sys.stdout.close()
                finally:
                    try:
                        sys.stderr.flush()
                    finally:
                        sys.stderr.close()
    return wrapper


@suppress_broken_pipe_msg
def main(source: ty.Union[str, ty.Iterable],
         out_fn: ty.Union[str, None],
         parser_options: dict,
         auto_config=False,
         skip_rows=None,
         skip_errors=True,
         max_errors=100,
         make_tabix: bool = False):
    try:
        parser = parsers.GenericGwasLineParser(**parser_options)
    except exceptions.ConfigurationException:
        parser = None
    finally:
        # Sniffer does not allow exact parser AND partial options to be specified. Once one is created, null the other.
        parser_options = None

    if source is None:
        source = sys.stdin

    if not auto_config and (skip_rows is None or parser is None):
        logger.error('Please provide all options required to parse the file, or use the --auto flag to guess')
        sys.exit(1)

    # Guess how to read the file. If no parser was provided, try to guess columns.
    reader = sniffers.guess_gwas_generic(source, skip_rows=skip_rows, parser=parser, parser_options=parser_options,
                                         skip_errors=skip_errors, max_errors=max_errors)

    try:
        dest_fn = reader.write(out_fn, make_tabix=make_tabix) or 'console'
    except exceptions.TooManyBadLinesException:
        logger.error('ERROR: Too many lines failed to parse; stopping.')
    except Exception:
        logger.exception('Conversion failed due to unknown error')
    else:
        logger.info('Conversion succeeded! Results written to: {}'.format(dest_fn))
    finally:
        for n, reason, _ in reader.errors:
            logger.error('Excluded row {} from output due to parse error: {}'.format(n, reason))


def run_cli():
    """Command line arguments"""
    parser = argparse.ArgumentParser(description="A utility to convert GWAS files into a standardized format."
                                                 "All column numbers should specify human-friendly indices (1-based)")
    parser.add_argument('file', nargs='?', default=None,
                        help='The path to the input GWAS file')
    parser.add_argument('--dest', type=str,
                        help='The filename to use when saving the output GWAS file')
    parser.add_argument('--skip-rows', dest='skip_rows', type=int, default=None,
                        help='The number of non-data lines that will be skipped')
    parser.add_argument('--max-errors', type=int, default=100, dest='max_errors',
                        help='How many unparseable lines can be skipped before the script stops')
    parser.add_argument('--stop-on-error', action='store_true', dest='stop_on_error',
                        help='Should we stop parsing immediately on the first error?')

    # Format specification options  # TODO: Improve argparse to make clear which options are mutually exclusive.
    # There are multiple ways of specifying chrom/pos/ref/alt columns: auto-detect, marker, or explicit columns
    parser.add_argument('--auto', action='store_true', dest='auto',
                        help='Automatically detect appropriate parser options')

    parser.add_argument('-m', '--marker_col', type=int, dest='marker_col',
                        help='Column number with marker (used to find chr, pos, ref, and alt)')

    parser.add_argument('-c', '--chrom_col', type=int, dest='chrom_col',
                        help='Column for chromosome')
    parser.add_argument('-b', '--pos_col', type=int, dest='pos_col',
                        help='Column for position')
    parser.add_argument('-r', '--ref_col', type=int, dest='ref_col',
                        help='Column for refvar')
    parser.add_argument('-a', '--alt_col', type=int, dest='alt_col',
                        help='Column for altvar')

    # This argument is always required
    parser.add_argument('-p', '--pvalue_col', type=int, dest='pvalue_col',
                        help='Column for pvalue (or logpvalue)')
    parser.add_argument('--is_neg_log_pvalue', action='store_true', dest='is_neg_log_pvalue',
                        help='Is the pvalue data stored in -log10(p) form?')

    # Optional fields
    parser.add_argument('--rsid_col', type=int, dest='rsid_col',
                        help='Column for rsid')
    parser.add_argument('--beta_col', type=int, dest='beta_col',
                        help='Column number for beta (effect size)')
    parser.add_argument('--stderr_beta_col', type=int, dest='stderr_beta_col',
                        help='Column for stderr of effect size')

    parser.add_argument('--allele_freq_col', type=int, dest='allele_freq_col',
                        help='Column for allele frequency')
    parser.add_argument('--allele_count_col', type=int, dest='allele_count_col',
                        help='Column for count of a specific allele (must also provide --n_samples_col)')
    parser.add_argument('--n_samples_col', type=int, dest='n_samples_col',
                        help='Column for number of samples (must also provide --allele_count_col')
    parser.add_argument('--is_ref_effect', action='store_true', dest='is_ref_effect',
                        help='By default, frequency will apply to alt allele. Specify this flag if it refers to ref.')

    parser.add_argument('--index', action='store_true', dest='index',
                        help='Whether to tabix and compress the output')

    args = parser.parse_args()
    parser_options = {
        fieldname: getattr(args, fieldname)
        for fieldname in ['marker_col', 'chrom_col', 'pos_col', 'ref_col', 'alt_col', 'pvalue_col',
                          'rsid_col', 'beta_col', 'stderr_beta_col',
                          'is_neg_log_pvalue', 'allele_freq_col', 'allele_count_col', 'n_samples_col']
    }
    parser_options['is_alt_effect'] = not args.is_ref_effect

    main(
        args.file,
        args.dest or None,
        parser_options,
        auto_config=args.auto,
        max_errors=args.max_errors,
        skip_errors=not args.stop_on_error,
        skip_rows=args.skip_rows,
        make_tabix=args.index
    )


if __name__ == '__main__':
    import time

    t1 = time.time()
    run_cli()
    logger.info('Conversion complete: Elapsed time: {}s'.format(time.time() - t1))
