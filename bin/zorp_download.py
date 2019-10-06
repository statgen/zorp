"""
Downloader script. Finds all assets available for download.

We can draw inspiration from the NLTK downloader script, here:
    https://github.com/nltk/nltk/blob/develop/nltk/downloader.py

The repository is defined as a static asset url with a JSON file, as follows:


{
  grch37: {
    genes: {
      version: '_',
      files: [
        { path: _, size: _, sha256: _ }
      ]
    },
    rsids: {
      version: '()',
      files: []
    },
  },
  grch38: {
    genes: ...
    rsids: ...
  }
}
"""
import argparse
import os
import urllib

DEFAULT_URL = 'https://example.dev'  # TODO: Change this to a CSG server!

# TODO: Make this setting consistent

os.environ.get()

def run_cli():
    parser = argparse.ArgumentParser(description='Download data assets for lookup features')
    parser.add_argument('--url', help='The URL of the download site specifying article information')
    pass


if __name__ == "__main__":
    pass
