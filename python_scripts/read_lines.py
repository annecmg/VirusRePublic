#!/usr/bin/env python3
"""Author: Devin van Valkengoed

Date: Sep-2024
Description: This module contains a function that can read lines from a .txt
             file and is imported in several Snakemake scripts in this repo.
"""
from snakemake.io import expand


def read_lines(input_file):
    """Returns a list of accession obtained from a given .txt file

    Key arguments:
        input_file -- str, full pathway to a file that contains accessions
                      on every new line.

    Returns:
        samples -- lst, containing accessions as elements.

    Note: the function will skipp empty lines and comment lines starting with
    '#'.
    """
    samples = []
    with open(input_file) as f:
        for sample in f.readlines():
            sample = sample.strip()
            if not sample:
                continue
            if sample.startswith("#"):
                continue
            if len(sample) > 0:
                samples.append(sample)
    f.close()

    return expand("{sample}", sample=samples)


if __name__ == '__main__':
    """This module is imported as a line reading helper module"""
    pass
