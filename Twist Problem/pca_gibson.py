# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 20:51:09 2023

@author: ji619
"""

import itertools

from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO


# Function to find duplicate regions
def find_duplicates(sequence1, sequence2, min_duplicate_length):
    duplicates = []
    for i in range(len(sequence1)):
        for j in range(i + min_duplicate_length, len(sequence1) + 1):
            subsequence1 = sequence1[i:j]
            if subsequence1 in sequence2:
                alignment = pairwise2.align.localms(subsequence1, sequence2, 2, -1, -1, -1)
                best_alignment = max(alignment, key=lambda x: x[2])
                if best_alignment[2] >= min_duplicate_length:
                    duplicates.append(subsequence1)
    return list(set(duplicates))


def PCA2(seqs, min_duplicate_length, max_recur=2):
    res = []

    def _PCA(seq_list, recur=1):
        if recur == max_recur:
            res.extend(seq_list)
            print(len(res), len(seq_list))
            return
        for s1, s2 in itertools.product(seq_list, seqs):
            _PCA(PCA(s1, s2, min_duplicate_length), recur=recur + 1)

    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if i < j:
                _PCA(PCA(seqs[i], seqs[j], min_duplicate_length))

    return res


# Function to process two DNA sequences and find duplicates
def PCA(sequence1, sequence2, min_duplicate_length):
    # Read sequences from files
    sequence2 = sequence2.reverse_complement()

    # Find duplicate regions
    duplicates = find_duplicates(sequence1, sequence2, min_duplicate_length)

    # Create sequence1a and sequence2a lists
    sequence1a_list = []
    sequence2a_list = []

    # Create sequence1a by ending in each duplicate value in sequence1
    for duplicate_value in duplicates:
        if duplicate_value in sequence1:
            end_index = sequence1.rindex(duplicate_value) + len(duplicate_value)
            sequence1a = sequence1[:end_index]
            sequence1a_list.append(sequence1a)

    # Create sequence2a by starting from each duplicate value in sequence2
    for duplicate_value in duplicates:
        if duplicate_value in sequence2:
            start_index = sequence2.index(duplicate_value) + len(duplicate_value)
            sequence2a = sequence2[start_index:]
            sequence2a_list.append(sequence2a)

    # Combine sequence1a and sequence2a
    combine = [a + b for a, b in zip(sequence1a_list, sequence2a_list)]

    return combine


def Gibson(sequence1, sequence2, min_duplicate_length):
    # Find duplicate regions
    duplicates = find_duplicates(sequence1, sequence2, min_duplicate_length)

    # Create sequence1a and sequence2a lists
    sequence1a_list = []
    sequence2a_list = []

    # Create sequence1a by ending in each duplicate value in sequence1
    for duplicate_value in duplicates:
        if duplicate_value in sequence1:
            end_index = sequence1.rindex(duplicate_value) + len(duplicate_value)
            sequence1a = sequence1[:end_index]
            sequence1a_list.append(sequence1a)

    # Create sequence2a by starting from each duplicate value in sequence2
    for duplicate_value in duplicates:
        if duplicate_value in sequence2:
            start_index = sequence2.index(duplicate_value) + len(duplicate_value)
            sequence2a = sequence2[start_index:]
            sequence2a_list.append(sequence2a)

    # Combine sequence1a and sequence2a
    combine = [a + b for a, b in zip(sequence1a_list, sequence2a_list)]

    return combine


if __name__ == '__main__':
    # Example usage PCA:
    sequence1 = Seq('TCCCTGGGCTCTTTTAGTGGACGGAGACCCAGCTGTCAGTTTGTTGTAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
    sequence2 = Seq('CTGCCCAAGCCTACCGTGAATCATCTAATCCCTCCATGGAGTAAGTGGTGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT')
    min_duplicate_length = 20  # Standard
    with open('data/sample3.fasta') as f:
        seqs = [s.seq for s in SeqIO.FastaIO.FastaIterator(f)][1:]

    combined_sequences = PCA2(seqs, min_duplicate_length)
    print(combined_sequences)

    # Example usage Gibson:
    min_duplicate_length = 15  # Standard

    combined_sequences = Gibson(sequence1, sequence2, min_duplicate_length)
    # print(combined_sequences)
