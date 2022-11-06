#!/usr/bin/env python3

#  Copyright (c) 2022. Dennis Wiersma.
#  Licensed under GPLv3. See LICENSE file.

"""
Module description
"""

# METADATA #

# IMPORTS #
import sys
import os
import argparse
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# CODE #
class MSAHandler:

    def __init__(self, path_to_protein_family):
        self.path_to_protein_family = path_to_protein_family

    def msa_from_fasta(self):
        try:
            msa_file = os.popen("/Users/denniswiersma/clustalo -i " + self.path_to_protein_family)
            return msa_file
        except ChildProcessError:
            print("Clustalo does not function correctly. Please check your install.")

    def read_msa(self):
        """
        Read the MSA file using BioPython's AlignIO object.

        :return: AlignIO object containing the MSA data.
        """
        try:
            return AlignIO.read(self.msa_from_fasta(), "fasta")
        except FileNotFoundError:
            print("Could not read MSA.")


def get_arguments():
    """
    Creates arguments and handles interfacing with the cli using an ArgumentParser object.

    :return: parsed arguments from the cli.
    """
    parser = argparse.ArgumentParser(prog="Smoking Gun",
                                     description="Determine SNP severity for a given protein family.")
    parser.add_argument("ProteinFile",
                        type=str,
                        help="Path to the fasta file containing the protein family.")
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

    msa_object = MSAHandler(args.ProteinFile)
    align_object = msa_object.read_msa()
    for record in align_object:
        print(record)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
