import time
import os
import motif_finder as finder
from scipy.special import rel_entr
import shutil
import numpy as np

result_relative_entropy = {}
result_overlapping_sites = {}
result_percent_overlap_sites = {}
result_overlapping_positions = {}
result_percent_overlap_positions = {}
result_runtime = {}

average_relative_entropy = {}
average_overlapping_sites = {}
average_percent_overlap_sites = {}
average_overlapping_positions = {}
average_percent_overlap_positions = {}
average_runtime = {}

# helper funtion to calculate the relative entropy between motif.txt and predictedmotif.txt
def Relative_Entropy(dataset_description):
    motif = import_motif('benchmarks/' + str(dataset_description) + '/motif.txt')
    predictedmotif = import_motif('results/' + str(dataset_description) + '/predictedmotif.txt')
    rel_ent = 0
    # print(motif)
    # print(predictedmotif)
    for i in range(len(motif)):
        # print(motif[i])
        # print(predictedmotif[i])
        row_diff = rel_entr(motif[i], predictedmotif[i])
        rel_ent += row_diff
    return sum(rel_ent)

#  helper function to import the motif file from the path
def import_motif(dataset_description):
	with open(dataset_description,'r') as f:
		lines = f.readlines()
	motif = []
	for line in lines[1:-1]:
		motif.append(list(map(float,line.strip().split())))
	return motif

# helper function to calculate the overlapping sites between experiment results and ground truth
# IMPORTANT!!!: Note that two sites overlap if at least ML/2 of their positions are common
def number_overlapping_sites(dataset_description):
    with open('benchmarks/' + str(dataset_description) + '/sites.txt', 'r') as f:
        sites = list(map(int, f.readlines()))
    with open('results/' + str(dataset_description) + '/predictedsites.txt', 'r') as f:
        predictedsites = list(map(int, f.readlines()))
    with open('benchmarks/' + str(dataset_description) + '/motiflength.txt', 'r') as f:
        motif_length = int(f.readlines()[0])

    overlap_sites = 0
    for i in range(len(sites)):
        gap = np.absolute(sites[i] - predictedsites[i])
        if gap <= motif_length // 2:
            overlap_sites += 1

    return overlap_sites, overlap_sites/len(sites)

# helper function to calculate the overlapping positions between experiment results and ground truth
def number_overlapping_positions(dataset_description):
    with open('benchmarks/' + str(dataset_description) + '/sites.txt', 'r') as f:
        sites = list(map(int, f.readlines()))
    with open('results/' + str(dataset_description) + '/predictedsites.txt', 'r') as f:
        predictedsites = list(map(int, f.readlines()))
    with open('benchmarks/' + str(dataset_description) + '/motiflength.txt', 'r') as f:
        motif_length = int(f.readlines()[0])

    overlap_positions = 0
    for i in range(len(sites)):
        gap = np.absolute(sites[i] - predictedsites[i])
        if gap <= motif_length:
            overlap_positions += motif_length - gap

    return overlap_positions, overlap_positions/(len(sites)*motif_length)



directory = "results"
parent_dir = os.getcwd()
path = os.path.join(parent_dir, directory)
if os.path.exists(path) and os.path.isdir(path):
    shutil.rmtree(path)
os.mkdir(path)

directory = "benchmarks"
parent_dir = os.getcwd()
path = os.path.join(parent_dir, directory)

#  These are list containing the average performance over the 10 data sets for each value of the hyper-parameters
relative_entropy = []
num_overlap_sites = []
percent_overlap_sites = []
num_overlap_positions = []
percent_overlap_positions = []
running_time = []
index = 0

for subdir, dirs, files in os.walk(path):
    for dir_name in dirs:

        # running time
        start_time = time.monotonic()
        finder.find_motif(dir_name)
        end_time = time.monotonic()
        delta = round((end_time - start_time), 3)
        result_runtime[dir_name] = delta
        print('{} completed! Running time: {}'.format(str(dir_name), str(delta)))

        # relative entropy (i.e.,KL divergence)
        rel_entropy = Relative_Entropy(dir_name)
        result_relative_entropy[dir_name] = rel_entropy

        # number amd percentage of overlapping sites
        number_sites, percent_sites = number_overlapping_sites(dir_name)
        result_overlapping_sites[dir_name] = number_sites
        result_percent_overlap_sites[dir_name] = percent_sites

        # number amd percentage of overlapping positions
        number_positions, percent_positions = number_overlapping_positions(dir_name)
        result_overlapping_positions[dir_name] = number_positions
        result_percent_overlap_positions[dir_name] = percent_positions

        running_time.append(delta)
        num_overlap_sites.append(number_sites)
        percent_overlap_sites.append(percent_sites)
        num_overlap_positions.append(number_positions)
        percent_overlap_positions.append(percent_positions)
        relative_entropy.append(rel_entropy)

        index += 1

        # get the average results once a set of 10 datasets are done
        if index == 10:
            average_relative_entropy[dir_name[0:-5]] = np.mean(relative_entropy)
            average_overlapping_sites[dir_name[0:-5]] = np.mean(num_overlap_sites)
            average_percent_overlap_sites[dir_name[0:-5]] = np.mean(percent_overlap_sites)
            average_overlapping_positions[dir_name[0:-5]] = np.mean(num_overlap_positions)
            average_percent_overlap_positions[dir_name[0:-5]] = np.mean(percent_overlap_positions)
            average_runtime[dir_name[0:-5]] = np.mean(running_time)
            relative_entropy = []
            num_overlap_sites = []
            percent_overlap_sites = []
            num_overlap_positions = []
            percent_overlap_positions = []
            running_time = []
            index = 0
            print("----------------------This set of 10 datasets is completed----------------------")



# print(result_relative_entropy)
# print(result_overlapping_sites)
# print(result_percent_overlap_sites)
# print(result_overlapping_positions)
# print(result_percent_overlap_positions)
# print(result_runtime)

# print(average_relative_entropy)
# print(average_overlapping_sites)
# print(average_percent_overlap_sites)
# print(average_overlapping_positions)
# print(average_percent_overlap_positions)
# print(average_runtime)
