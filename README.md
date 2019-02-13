# ZORP: A helpful GWAS parser

## Usage
```python
import readers
import parsers

# Create a reader instance
reader = readers.TabixReader('input.bgz', parser=parsers.standard_gwas_parser)
reader.add_filter('chr', '19')
reader.add_filter('pvalue', lambda val, line: val < 5e-8)

for row in reader.fetch('X', 1, 1_000_000):
    print(row)
    
for row in reader:
    print(row.chrom)
```


## Development

To run unit tests, use
`PYTHONPATH=. pytest tests --flake8`






