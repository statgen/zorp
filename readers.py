"""
External-facing API to load GWAS data from files
"""


def read_gwas(filename:str, **kwargs):
    """Read a GWAS file, guessing the dialect if necessary"""
    # TODO: Find a mechanism where any new subclasses could automatically be used by the dialect detection mechanism?
    pass




# TODO Future utils brainstorming
#  -- How to handle if a gwas is divided into files? (maybe a sniffer with a directory-detection mode that chains iterators and verifies chain order?)
#  -- Alignment and processing logic (across files? What would we want the output to be then?)
# Experiment with parsing more files. Can we do this async/ multiprocessing? What would that get us? eg coroutines for the iteration and parsing steps?
# Performance timings? Various approaches?