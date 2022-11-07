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
from Bio.Seq import MutableSeq


# CODE #

def align_from_fasta(path_to_protein_family):
    """
    Use Clustalo on the command line to open a fasta file and align it.

    :return: MSA data in an AlignIO object.
    """
    try:
        msa_data = os.popen("/Users/denniswiersma/clustalo -i " + path_to_protein_family)
    except ChildProcessError:
        print("Clustalo does not function correctly. Please check your install.")
        sys.exit()

    try:
        return AlignIO.read(msa_data, "fasta")
    except FileNotFoundError:
        print("Could not read MSA.")
        sys.exit()


def define_scoring(msa_data, snp_position):
    """
    Calculates the AA frequency at the position at which the SNP will be placed.
    Then normalises these values by dividing by the sum of AA frequencies at that position.
    :param msa_data: Data from the MSA as an AlignIO object
    :param snp_position: Position at which the SNP should be placed
    :return: A dictionary containing normalised frequency scores
    """
    # Initialise empty dict
    aa_frequency = {}

    # Loop through record of AlignIO object
    for record in msa_data:
        # Check if AA is in the dict already
        if record.seq[snp_position] in aa_frequency:
            # Add one to the count for that AA
            aa_frequency[record.seq[snp_position]] += 1
        else:
            # Set the count for that AA to one
            aa_frequency[record.seq[snp_position]] = 1

    # Calculate how many AA were counted in the entire dictionary
    sum_of_frequencies = sum(aa_frequency.values())

    # Normalise counts by dividing every dictionary value by the sum of frequencies calculated above
    for frequency in aa_frequency:
        aa_frequency[frequency] = float(aa_frequency[frequency]/sum_of_frequencies)

    return aa_frequency


def read_dna(dna_seq_or_file):
    """
    Takes either a DNA sequence directly or a fasta file and extracts it.

    :return: A SecRecord object.
    """
    # Check if sequence or filename was entered
    if all(i.upper() in "ATCG" for i in dna_seq_or_file):

        # Prompt user for said parameters
        print("Looks like you've entered a DNA sequence. "
              "You can optionally provide some more data. You may leave these empty if you wish.")
        seq_id = input("Please enter a sequence id: ")
        seq_name = input("Please enter a sequence name: ")
        seq_desc = input("Please enter a sequence description: ")

        # Return a SeqRecord instance
        return SeqRecord(
            MutableSeq(dna_seq_or_file),
            id=seq_id,
            name=seq_name,
            description=seq_desc
        )

    else:
        try:
            # Open fasta file containing DNA sequence
            with open(dna_seq_or_file) as dna_file:
                seq_record_immutable = list(SeqIO.parse(dna_file, "fasta"))[0]
                return SeqRecord(
                        MutableSeq(seq_record_immutable.seq),
                        id=seq_record_immutable.id,
                        name=seq_record_immutable.name,
                        description=seq_record_immutable.description
                )
        except FileNotFoundError:
            print("Could not read DNA sequence file.")
            sys.exit()


def insert_snp(protein_dna_record, snp_position, snp_letter):
    """
    Inserts an SNP into a DNA sequence extracted from a SeqRecord with MutableSeq.

    :param protein_dna_record: SeqRecord object containing the DNA to insert the SNP into
    :param snp_position: Position at which the SNP should be placed
    :param snp_letter: Letter that should replace the letter already present
    :return:
    """
    protein_dna = protein_dna_record.seq

    if snp_position == 0:
        protein_dna = snp_letter + protein_dna

    elif snp_position <= len(protein_dna):
        protein_dna[snp_position - 1] = snp_letter

    elif snp_position > len(protein_dna):
        protein_dna += snp_letter

    else:
        raise IndexError("Could not resolve SNP position.")

    return SeqRecord(
        MutableSeq(protein_dna),
        id=protein_dna_record.id,
        name=protein_dna_record.name,
        description=protein_dna_record.description
    )


def translate_dna_to_aa(dna_record):
    """
    Translates dna to amino acid using SeqRecord objects
    :param dna_record: SeqRecord object containing the DNA sequence
    :return:
    """
    protein_seq = dna_record.translate()
    return SeqRecord(
        MutableSeq(protein_seq.seq),
        id=dna_record.id,
        name=dna_record.name,
        description=dna_record.description
    )


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
                        nargs=2,
                        help="Location of the SNP in the DNA sequence and the nucleotide to change to.\n"
                             "Example: 6 A\n"
                             "This would change position 6 to an A.")
    parser.add_argument("-s", "--sequence",
                        required=True,
                        type=str,
                        help="DNA sequence containing the SNP."
                             "Can be either the path to a FASTA file or a DNA sequence directly.")
    return parser.parse_args()


def main(args):
    """Function description"""
    # Get arguments from the cli
    arguments = get_arguments()

    # Fetch fasta data and perform msa
    msa_data = align_from_fasta(arguments.ProteinFile)

    # Fetch DNA sequence
    protein_dna = read_dna(arguments.sequence)

    # Insert SNP into DNA sequence
    snp_position = int(arguments.snp[0])
    snp_letter = arguments.snp[1]
    protein_dna_snp = insert_snp(protein_dna, snp_position, snp_letter)

    # Translate DNA with SNP to protein
    protein_sequence_snp = translate_dna_to_aa(protein_dna_snp)

    # Generate scoring based on the msa
    msa_based_scoring = define_scoring(msa_data, snp_position)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
