"""
This is the program to generate the benchmarks with following hyper-parameters

Let the default parameter combination be “ICPC = 2, ML = 8, SL = 500, SC = 10”
(a) set ICPC = 1, 1.5, 2, while other parameters are at default values.
(b)) Set ML = 6,7,8, while other parameters are at default values.
(c) Set SC = 5, 10, 20, while other parameters are at default values.
You may keep SL fixed at 500.
"""

import numpy as np
import os
import random
import shutil

default_ML = 8
default_SC = 10
default_ICPC = 2
default_SL = 500
nucleotides = ['A', 'C', 'G', 'T']

ICPC = [1.0, 1.5, 2.0]
ML = [6, 7]
SC = [5, 20]

# helper function to generate SC random sequences (with uniform nucleotide frequencies).
# Each random sequence has length SL.
def generate_sequence(sequence_length):
    out=''
    for i in range(sequence_length):
        out += nucleotides[int(random.randint(0,100000)%4)]
    return out

# helper function to generate a random motif (position weight matrix) of length ML
# with total information content (as discussed in class) being ICPC * ML
# As mentioned in the doc:
#     p ~ 0.8105 ( ICPC ~ 1.0)
#     p ~ 0.9245 ( ICPC ~ 1.5)
#     p = 1.0 ( ICPC = 2.0)
def generate_motif(motif_length, ICPC):
    motif = []
    if ICPC == 2:
        for col in range(motif_length):
            preferred_nucleotide = int(random.randint(0,100000)%4)
            column = [0, 0, 0, 0]
            column[preferred_nucleotide] = 1
            motif.append(column)
    else:
        for col in range(motif_length):
            if ICPC == 1.5:
                prob_preferred_nucleotide = 0.9245
            else:
                prob_preferred_nucleotide = 0.8105
            prob_others = (1 - prob_preferred_nucleotide) / 3
            preferred_nucleotide = int(random.randint(0, 100000) % 4)
            column = [0, 0, 0, 0]
            for i in range(4):
                if i == preferred_nucleotide:
                    column[i] = prob_preferred_nucleotide
                else:
                    column[i] = prob_others
            motif.append(column)

    return np.array(motif)

# helper function to generate a pattern by sampling from this random motif
def generate_motif_pattern(motif_length, pwm):
    motif_pattern = ''
    for i in range(motif_length):
        motif_pattern += np.random.choice(nucleotides, p=pwm[i])
    return motif_pattern

# helper function to “Plant” one sampled site at a random location in each random sequence
def plant_motif(motif_length, motif_pattern, sequence):
    site = np.random.randint(0, len(sequence) - motif_length)
    end_position = site + motif_length
    string = sequence[0:site] + motif_pattern + sequence[end_position:]
    return string, site

# helper function to write to the files
def write_to_file(dataset_description, sequences, sites, pwm, motif_length):
    os.mkdir('benchmarks/' + str(dataset_description))
    f = open('benchmarks/' + str(dataset_description) + '/sequences.fa', 'w+')
    for index, sequence in enumerate(sequences):
        f.write('>sequence' + str(index+1) + '\n')
        f.write(sequence + '\n')
    f.close()
    f = open('benchmarks/' + str(dataset_description) + '/sites.txt', 'w+')
    for site in sites:
        f.write(str(site) + '\n')
    f.close()
    f = open('benchmarks/' + str(dataset_description) + '/motif.txt', 'w+')
    f.write('>MOTIF {}\t {}\n'.format(str(dataset_description), str(motif_length)))
    for col in pwm:
        column = ' '.join(map(str, col))
        f.write(column + '\n')
    f.write('<')
    f.close()
    f = open('benchmarks/' + str(dataset_description) + '/motiflength.txt', 'w+')
    f.write(str(motif_length))
    f.close()


directory = "benchmarks"
parent_dir = os.getcwd()
path = os.path.join(parent_dir, directory)
if os.path.exists(path) and os.path.isdir(path):
    shutil.rmtree(path)
os.mkdir(path)

# generates 10 data sets for each hyper-parameters setting
# There are (1+2+2+2) x 10 = 70 data sets in the benchmark
for a in range(10):
    for x in ICPC:
        sequences = []
        motif_patterns = []
        sites = []
        pwm = generate_motif(default_ML, x)
        for i in range(default_SC):
            motif_pattern = generate_motif_pattern(default_ML, pwm)
            sequence = generate_sequence(default_SL)
            modified_seq, site = plant_motif(default_ML, motif_pattern, sequence)
            motif_patterns.append(motif_pattern)
            sequences.append(modified_seq)
            sites.append(site)

        dataset_description = ("ICPC_" + str(x) + "_NO_" + str(a))
        write_to_file(dataset_description, sequences, sites, pwm, default_ML)

    for motif_length in ML:
        sequences = []
        motif_patterns = []
        sites = []
        pwm = generate_motif(motif_length, default_ICPC)
        for i in range(default_SC):
            motif_pattern = generate_motif_pattern(motif_length, pwm)
            sequence = generate_sequence(default_SL)
            modified_seq, site = plant_motif(motif_length, motif_pattern, sequence)
            motif_patterns.append(motif_pattern)
            sequences.append(modified_seq)
            sites.append(site)

        dataset_description = ("ML_" + str(motif_length) + "_NO_" + str(a))
        write_to_file(dataset_description, sequences, sites, pwm, motif_length)

    for sequence_count in SC:
        sequences = []
        motif_patterns = []
        sites = []
        pwm = generate_motif(default_ML, default_ICPC)
        for i in range(sequence_count):
            motif_pattern = generate_motif_pattern(default_ML, pwm)
            sequence = generate_sequence(default_SL)
            modified_seq, site = plant_motif(default_ML, motif_pattern, sequence)
            motif_patterns.append(motif_pattern)
            sequences.append(modified_seq)
            sites.append(site)

        dataset_description = ("SC_" + str(sequence_count) + "_NO_" + str(a))
        write_to_file(dataset_description, sequences, sites, pwm, default_ML)

