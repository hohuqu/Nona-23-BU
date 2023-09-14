# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 20:51:09 2023

@author: ji619
"""

from Bio import pairwise2

# Function to find the complement of a DNA sequence
def find_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement_dict[base] for base in sequence)

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

# Function to read a DNA sequence from a file
def read_sequence_from_file(filename):
    with open(filename, 'r') as file:
        return file.read().strip()

# Function to process two DNA sequences and find duplicates
def PCA(file1, file2, min_duplicate_length):
    # Read sequences from files
    sequence1 = read_sequence_from_file(file1)
    sequence2 = find_complement(read_sequence_from_file(file2))[::-1]

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

# Example usage PCA:
file1 = "oligo.txt"
file2 = "oligo2.txt"
min_duplicate_length = 20 #Standard

combined_sequences = PCA(file1, file2, min_duplicate_length)
print(combined_sequences)

def Gibson(file1, file2, min_duplicate_length):
    # Read sequences from files
    sequence1 = read_sequence_from_file(file1)
    sequence2 = read_sequence_from_file(file2)

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

# Example usage Gibson:
file1 = "oligo.txt"
file2 = "oligo2.txt"
min_duplicate_length = 15 #Standard

combined_sequences = Gibson(file1, file2, min_duplicate_length)
print(combined_sequences)
