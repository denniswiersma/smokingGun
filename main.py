#!/usr/bin/env python3

#  Copyright (c) 2022. Dennis Wiersma.
#  Licensed under GPLv3. See LICENSE file.

"""
Module description
"""

# METADATA #

# IMPORTS #
import sys
from Bio import AlignIO


# CODE #

def main(args):
    """Function description"""
    alignment = AlignIO.read(open("galectin3/msa.fasta"), "fasta")
    for record in alignment:
        print(record)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
