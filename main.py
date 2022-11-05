#!/usr/bin/env python3

#  Copyright (c) 2022. Dennis Wiersma.
#  Licensed under GPLv3. See LICENSE file.

"""
Module description
"""

# METADATA #

# IMPORTS #
import sys
import argparse
from Bio import AlignIO


# CODE #
class MSAHandler:

    def __init__(self, path_to_msa):
        self.path_to_msa = path_to_msa

    def read_msa(self, path_to_msa):
        """
        Read the MSA file using BioPython's AlignIO object.

        :param path_to_msa: filepath to the file containing the MSA.
        :return: AlignIO object containing the MSA data.
        """
        try:
            with open(path_to_msa) as msa_file:
                return AlignIO.read(msa_file, "fasta")
        except FileNotFoundError:
            print("Could not find file. Please try again.")


def get_arguments():
    """
    Creates arguments and handles interfacing with the cli using an ArgumentParser object.

    :return: parsed arguments from the cli.
    """
    parser = argparse.ArgumentParser(prog="Smoking Gun",
                                     description="Determine SNP severity for a given MSA.")
    parser.add_argument("MSAfile",
                        type=str,
                        help="Path to the MSA.")
    parser.add_argument("-p", "--snp",
                        required=True,
                        type=int,
                        help="Location of the SNP in the DNA sequence.")
    parser.add_argument("-s", "--sequence",
                        required=True,
                        type=str,
                        help="DNA sequence containing the SNP."
                             "Can be either the path to a FASTA file or a DNA sequence directly.")
    return parser.parse_args()


def main(args):
    """Function description"""
    args = get_arguments()

    msa_object = MSAHandler(args.MSAfile)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
