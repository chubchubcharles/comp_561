import matplotlib.pyplot as plt
import math
import numpy as np
import random
import scipy
import pylab
import csv
from sets import Set
import scipy.cluster.hierarchy as sch
from Bio.Cluster import *
from Bio import Phylo
from Bio.Cluster import Node, Tree
from Bio.Cluster import distancematrix,treecluster
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list

#var declarations
datasets = []
dist_mat = np.zeros((3,3))
position = ""
chromosome = 1

#methods
def create_plot():
	plt.subplot(3, 1, 1)
	plt.title('Random Copy Number Frequency Across Genome')
	create_subplots(True)
	plt.subplot(3, 1, 2)
	create_subplots(True)
	plt.subplot(3, 1, 3)
	create_subplots(False)
	plt.xlabel('Genomic Position')
	plt.show()

def create_subplots(pRandom):
	list = create_data(pRandom)
	datasets.append(list)
	#plot_data(list)

def create_data(pRandom):
	list = []
	for x in xrange(0, 20):
		if pRandom:
			list.append(random.uniform(-2.0,2.0))
		else:
			list.append(0)
	return list
	
def plot_data(pList):
	plt.plot(pList, color = 'black')
	plt.ylabel('Copy Number Frequency')
	plt.axis([-2, 22, -3, 3])

def fill_dist_mat():
	
	for row in xrange(0,len(dist_mat)):
		for col in xrange(0,row):
			dist_mat[row][col] = calc_sep_score(datasets[row], datasets[col])
			print (dist_mat[row][col])
	
	# dist_mat = symmetrize(dist_mat)

# def symmetrize(pMatrix):
# 	return pMatrix + pMatrix.T - np.diag(pMatrix.diagonal())

def calc_sep_score(vect1, vect2):
	squareSum =0
	size= len(vect1)
	for pos in xrange (0, size):
		squareSum += math.pow ((vect2[pos]- vect1[pos]),2.0)
	distance = math.sqrt(squareSum)
	return distance
	
def compute_dist_Mat (dataset, metric):
	dist = pdist(dataset, metric= metric)
	return dist

def print_matrix():
	for row in range(len(dist_mat)):
		for col in range(row):
			print dist_mat[row][col]

#main method
def main(numbofdata):
	# for pos in xrange (0,numbofdata):
	# 	create_subplots(True)
	#print datasets
	#create_subplots(True)
	#fill_dist_mat()
	#print dist_mat
	"""distance= distancematrix(datasets)
	tree = treecluster(distancematrix=distance)
	tree.scale()"""
	#plt.plot()
	
	distMat = compute_dist_Mat(datasets,'euclidean')
	Z = sch.linkage(distMat, method='average')
	leaves = leaves_list(Z)	
	numleaves= len(leaves)
	plt.subplot(1, 2, 1)
	leaf_labels = ['thyroid', 'cervix', 'large_intestine', 'stomach', 
	'central_nervous_system', 'haematopoietic_and_lymphoid_tissue', 'pancreas', 
 'urinary_tract', 'lung', 'breast', 'endometrium', 
 'skin', 'ovary', 'NS', 'prostate', 'kidney', 'liver']
	dendrogram(Z, orientation='right', labels = leaf_labels)
	for data in leaves:
		plt.subplot (numbofdata, 2, (2*data) +2)
		plt.plot (datasets[leaves_list(Z)[numleaves - data - 1]])
		# print data
		#plt.ylabel('log ratio number alteration')
	
	plt.suptitle('Dendogram from clustering of expression on chr' + str(chromosome))
	plt.show()

def parsingTSVFile ():

	cancer = set([])
	thyrCancer= []
	cervixCancer=[]
	LICancer=[] 
	StoCancer=[] 
	CNSCancer=[]
	HALCancer=[]
	pancreaCancer=[]
	UTCancer=[]
	lungCancer=[]
	breastCancer=[]
	endometriumCancer=[] 
	skinCancer=[]
	ovaryCancer=[]
	NSCancer=[]
	prostateCancer=[]
	kidneyCancer=[]
	liverCancer=[]
	#LintestineCancer stomachCancer CNSCancer haa
	dictCancer= {'thyroid':thyrCancer, 'cervix':cervixCancer, 'large_intestine':LICancer, 'stomach':StoCancer, 
	'central_nervous_system':CNSCancer, 'haematopoietic_and_lymphoid_tissue':HALCancer, 'pancreas':pancreaCancer, 
 'urinary_tract':UTCancer, 'lung':lungCancer, 'breast':breastCancer, 'endometrium':endometriumCancer, 
 'skin':skinCancer, 'ovary':ovaryCancer, 'NS':NSCancer, 'prostate':prostateCancer, 'kidney':kidneyCancer, 'liver':liverCancer}
	with open("CosmicCompleteCNA.tsv") as tsv:
		for line in csv.reader (tsv, dialect="excel-tab"):
			global chromosome 
			if line[12][:2]== str(chromosome) + ':':
				cancerType = line[3]
				cnv = line[8]
				#print cnv
				if cnv !='TOTAL_CN' and cnv!='':
					#print cnv
					try:
						sample_size = 12
						if len(dictCancer[cancerType]) < sample_size:
							cnv= float(cnv)
							dictCancer[cancerType].append(line)
					except ValueError as detail:
						print 'The error is due to:-----' , cnv, '---------',cnv == ''
					
					#cnv= float(cnv)
					#dictCancer[cancerType].append(cnv)
			#print line[12][0]
			cancer.add (line[3])
	# print len(cancer)
	# print dictCancer['cervix']

	#sorting every expression level in a given cancer
	#to parse for the start position: 
		# dictCancer['cervix'][1][12].split(':')[1].split("..")[0]
	for cancer in dictCancer:
		#sorting every expression level in a given cancer by genomic position
		dictCancer[cancer].sort(key = lambda x: x[12].split(':')[1].split("..")[0])
		#replacing dictionary values with expression levels only
		dictCancer[cancer] = [line[8] for line in dictCancer[cancer]]
		# appending to datasets
		datasets.append(dictCancer[cancer])

#levels(4)
parsingTSVFile ()
# print dictCancer['NS']
# print datasets
# print position
# limitations: select for chr 1, for 12 samples each, sorted by genomic location, 
main(len(datasets))


