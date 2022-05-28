################################################################################################

# Greedy Motif Search Algorithm to predict motif (position weight matrix) and predict sites in DNA sequences
# Input:  motiflength.txt
# 	      sequences.fa
# Output: predictedmotif.txt
#         predictedsites.txt

################################################################################################

import numpy as np
import os

motif_length = 0
sequences = []

# main method to find the common motif
def find_motif(dataset_description):
    with open('benchmarks/' + str(dataset_description) + '/motiflength.txt', 'r') as f:
        length = f.readlines()
    with open('benchmarks/' + str(dataset_description) + '/sequences.fa', 'r') as f:
        sequences = f.read().splitlines()
    sequences = [elem for elem in sequences if len(elem) == 500]
    motif_length = int(length[0])
    best_motif = np.zeros(len(sequences))
    positions = np.zeros(len(sequences))
    sequence_length = len(sequences[0])

    # Find two l-mers in sequences 1 and 2, form 2 x l  alignment matrix  and compute Score(s)
    # Pick the s with the highest score (i.e, information content)
    best_score = 0
    for i in range(sequence_length-motif_length + 1):
        positions[0] = i
        for j in range(sequence_length-motif_length + 1):
            positions[1] = j
            current_score = score(positions, 2, sequences, motif_length)
            # print(str(positions[0]) + " and " + str(positions[1]))
            # print(current_score)
            if current_score > best_score:
                best_motif[0] = positions[0]
                best_motif[1] = positions[1]
                best_score = current_score

    # print("---------------best-------------------")
    # print(str(best_motif[0]) + " and " + str(best_motif[1]))
    # print(best_score)

    positions[0] = best_motif[0]
    positions[1] = best_motif[1]

    # Iteratively add one l-mer from each of the other (t-2) sequences
    # At each of the following t-2 iterations, finds a “best” l-mer in sequence i
    for i in range(2, len(sequences)):
        best_score = 0
        for j in range(sequence_length-motif_length+1):
            positions[i] = j
            current_score = score(positions, i+1, sequences, motif_length)
            if current_score > best_score:
                best_score = current_score
                best_motif[i] = j
        positions[i] = best_motif[i]

    # print(positions)
    # print(best_motif)

    # “predictedsites.txt”, in the same format as “sites.txt” from the data set
    os.mkdir('results/' + str(dataset_description))
    f = open('results/' + str(dataset_description) + '/predictedsites.txt', 'w+')
    for i in range(len(best_motif)):
        f.write(str(int(best_motif[i])) + '\n')
    f.close()

    # “predictedmotif.txt”, in the same format as “motif.txt” from the data set
    f = open('results/' + str(dataset_description) + '/predictedmotif.txt', 'w+')
    f.write('>PMOTIF {}\t {}\n'.format(str(dataset_description), str(motif_length)))
    pwm = generate_pwm(best_motif, sequences, motif_length)
    for i in range(len(pwm)):
        f.write(str(pwm[i, 0]) + " " + str(pwm[i, 1]) + " " + str(pwm[i, 2]) + " " + str(pwm[i, 3]) + "\n")
    f.write("<")
    f.close()

# ----------------------------------------below are helper functions-----------------------------------------------

# helper function to give a score to the current pwm profile based on information content
def score(seq, num_seq, sequences, motif_length):
    profile_matrix = np.zeros((4, motif_length))
    information_content = 0
    for i in range(motif_length):
        for j in range(num_seq):
            if sequences[j][int(seq[j])+i] == 'A':
                profile_matrix[0, i] += 1
            elif sequences[j][int(seq[j])+i] == 'C':
                profile_matrix[1, i] += 1
            elif sequences[j][int(seq[j])+i] == 'G':
                profile_matrix[2, i] += 1
            elif sequences[j][int(seq[j])+i] == 'T':
                profile_matrix[3, i] += 1

    profile_matrix = profile_matrix/num_seq #Normalizes the profile matrix

    for k in range(motif_length):
        for x in range(4):
            if profile_matrix[x, k] != 0:
                information_content += profile_matrix[x, k] * np.log2(profile_matrix[x, k] / 0.25) # Information content function

    return information_content

# helper function to generate the pwm
def generate_pwm(best_motif, sequences, motif_length):
    profile_matrix = np.zeros((4, motif_length))
    for i in range(motif_length):
        for j in range(len(sequences)):
            if sequences[j][int(best_motif[j])+i] == 'A':
                profile_matrix[0, i] += 1
            elif sequences[j][int(best_motif[j])+i] == 'C':
                profile_matrix[1, i] += 1
            elif sequences[j][int(best_motif[j])+i] == 'G':
                profile_matrix[2, i] += 1
            elif sequences[j][int(best_motif[j])+i] == 'T':
                profile_matrix[3, i] += 1

    return profile_matrix.transpose()/len(sequences) # Outputs the normalized matrix

# find_motif('ICPC_2.0_NO_3')
