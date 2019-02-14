# ZORP: A helpful GWAS parser

## Why?
ZORP is intended to abstract away differences in file formats, and help you work with GWAS data from many 
different sources.

- Provide a single unified interface to read text, gzip, or tabixed data
- Separation of concerns between reading and parsing (with parsers that can handle the most common file formats)
- Includes helpers to auto-detect data format, filter for variants of interest, 

## Why not?
ZORP provides a high level abstraction. This means that it is convenient, at the expense of speed.

For GWAS files, ZORP does not sort the data for you, because doing so in python would be quite slow. You will still 
need to do some basic data preparation before using.

## Usage
```python
import readers
import parsers

# Create a reader instance
reader = readers.TabixReader('input.bgz', parser=parsers.standard_gwas_parser)

# We can filter data to the variants of interest. If you use a domain specific parser, columns can be referenced by name
reader.add_filter('chrom', '19')
reader.add_filter('pvalue', lambda val, line: val < 5e-8)

# Tabix files support iterating over all or part of the file
for row in reader.fetch('X', 500_000, 1_000_000):
    print(row)

# Iteration returns tuples (or namedtuples) of cleaned, parsed data
for row in reader:
    print(row.chrom)

# Write a compressed, tabix-indexed file containing the subset of variants of interest, in the parser's normalized format
out_fn = reader.write('outfile.txt', ['chrom', 'pos', 'pvalue'], make_tabix=True)
```


## Development

To run unit tests, use
`PYTHONPATH=. pytest tests --flake8`






