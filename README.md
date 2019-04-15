# ZORP: A helpful GWAS parser

## Why?
ZORP is intended to abstract away differences in file formats, and help you work with GWAS data from many 
different sources.

- Provide a single unified interface to read text, gzip, or tabixed data
- Separation of concerns between reading and parsing (with parsers that can handle the most common file formats)
- Includes helpers to auto-detect data format and filter for variants of interest 

## Why not?
ZORP provides a high level abstraction. This means that it is convenient, at the expense of speed.

For GWAS files, ZORP does not sort the data for you, because doing so in python would be quite slow. You will still 
need to do some basic data preparation before using.

## Usage
### Python
```python
from zorp import readers, parsers

# Create a reader instance
reader = readers.TabixReader('input.bgz', parser=parsers.standard_gwas_parser, skip_rows=1, skip_errors=True)

# We can filter data to the variants of interest. If you use a domain specific parser, columns can be referenced by name
reader.add_filter('chrom', '19')
reader.add_filter('log_pvalue', lambda val, line: val > 7.301)

# Iteration returns tuples (or namedtuples) of cleaned, parsed data
for row in reader:
    print(row.chrom)

# Tabix files support iterating over all or part of the file
for row in reader.fetch('X', 500_000, 1_000_000):
    print(row)

# Write a compressed, tabix-indexed file containing the subset of variants that match filters, choosing only specific 
#   columns. The data written out will be cleaned and standardized by the parser into a well-defined format. 
out_fn = reader.write('outfile.txt', columns=['chrom', 'pos', 'pvalue'], make_tabix=True)

# Real data is often messy. If a line fails to parse, the problem will be recorded.
for number, message, raw_line in reader.errors:
    print('Line {} failed to parse: {}'.format(number, message))

```

### Command line file conversion
The file conversion feature of zorp is also available as a command line utility. See `zorp-convert --help` for details.

This utility is currently in beta; please inspect the results carefully.

To auto-detect columns based on a library of commonly known file formats:

`$ zorp-convert --auto infile.txt --dest outfile.txt --compress`

Or specify your data columns exactly: 

`$ zorp-convert infile.txt --dest outfile.txt --compress  --skip-rows 1 --chr_col 1 --pos_col 2 --ref_col 3 --alt_col 4 pval_col 5`


## Development

To install dependencies and run in development mode:

`pip install -e '.[test]'`

To run unit tests, use

```bash
$ flake8 zorp
$ mypy zorp
$ pytest tests/
```
