import matplotlib.pyplot as plt
import numpy as np
import random
from Bio.Cluster import *
from Bio import Phylo
import pylab
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list

#var declarations
numData = 0
datasets = []
dist_mat = np.zeros((numData,numData))

#methods

def create_data(pNumSeq):
	global numData 
	numData = pNumSeq
	for x in range(pNumSeq-1):
		create_seq(True)
	create_seq(False)

def create_seq(pRandom):
	list = []
	for x in range(0, 20):
		if pRandom:
			list.append(random.uniform(-2.0,2.0))
		else:
			list.append(0)
	datasets.append(list)

def create_graph():
	if (len(datasets) != 0):
		for graphs in range(numData):
			plt.subplot(numData, 1, graphs + 1)
			plot_data(datasets[graphs])

		plt.ylabel('Expression Level')
		plt.xlabel('Genomic Position')
		plt.suptitle('Copy Number Frequency Across Genome', fontsize=12)
		plt.gcf().canvas.set_window_title('Copy Number Frequency Across Genome')
		plt.show()
	else :
		print "Datasets not created yet."

def plot_data(pList):
	plt.plot(pList, color = 'black')
	# plt.ylabel('Copy Number Frequency')
	plt.axis([-2, 22, -3, 3])

def calc_my_dist_mat():
	for row in range(len(dist_mat)):
		for col in range(row):
			dist_mat[row][col] = calc_my_dist_score(datasets[row], datasets[col])
	
	# dist_mat = symmetrize(dist_mat)

# def symmetrize(pMatrix):
# 	return pMatrix + pMatrix.T - np.diag(pMatrix.diagonal())

def calc_my_dist_score(pList1, pList2):
	sum = 0
	for index in range(len(pList1)):
		sum += (pList1[index] - pList2[index])
	return sum

def print_my_dist_mat():
	for row in range(len(dist_mat)):
		for col in range(row):
			print "difference between row: " + str(row) + " col: " + str(col) + " is " + str(dist_mat[row][col])

def calc_dist_mat():
	#distance matrix, calculated with Euclidean distance, even though default = 'e'
	matrix = distancematrix(datasets, dist='e')
	print "biopython distance matrix-e: " + str(matrix)

def create_dendrogram():
	plt.subplot(121)
	plt.ylabel('Samples')
	plt.xlabel('Distance')
	data_dist = pdist(datasets) # computing the distance
	data_link = linkage(data_dist, method='single') # computing the linkage
	dendrogram(data_link, orientation='right')

	for graph in range(numData):
		plt.subplot(numData, 2, (2*graph) + 2)
		plot_data(datasets[leaves_list(data_link)[numData - graph - 1]])

	# plt.subplot(322)
	# plot_data(datasets[leaves_list(data_link)[2]])
	# plt.subplot(324)
	# plot_data(datasets[leaves_list(data_link)[1]])
	# plt.subplot(326)
	# plt.xlabel('Genomic Position')
	# plot_data(datasets[leaves_list(data_link)[0]])
	# print str(leaves_list(data_link))

	plt.gcf().canvas.set_window_title('Hierarchical Clustering of Samples')
	plt.suptitle('Hierarchical Clustering of Samples', fontweight='bold', fontsize=14);
	plt.show()

#main method
create_data(5)
# create_graph()
create_dendrogram()
