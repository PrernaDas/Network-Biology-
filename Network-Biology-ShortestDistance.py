#!/usr/bin/python

import networkx as nx
import matplotlib.pyplot as plt
import random
import sys
from collections import deque
import timeit
import statistics
import math


start = timeit.default_timer()

# 1. Simplified Reactome Interaction Network

# nodes and edge list in the simplified Reactome Interaction Network
edge_list_ppisg = open("edge_list_ppisg.txt", "r")

nodes_ppisg = open("nodes_ppisg.txt", "r")


# Step-1: Build the Reactome Interaction Network graph which has many components
PPI = nx.Graph()
nodes = []
for n in nodes_ppisg:
	n = n.rstrip()
	nodes.append(n)
	
edges = []
for line in edge_list_ppisg:
	edges.append(line)

for edge in edges:
		edge = edge.split("\t")
		edge = edge[0], edge[1].rstrip()
		PPI.add_edge(edge[0], edge[1])
		
for n in nodes:
	if n not in PPI:
		PPI.add_node(n)

print "The Reactome Interaction Network, with many components", nx.info(PPI) # 6890 nodes		

#-------------------------------------------------------------------------------------------------------------#--------------#

# Step-2: Largest component of the Reactome Interaction Network	
PPI_largest_comp = max(nx.connected_component_subgraphs(PPI), key=len)
print "The largest component of the Reactome Interaction Network:", nx.info(PPI_largest_comp) # 6835 nodes	

#-------------------------------------------------------------------------------------------------------------#
# Step-3: Map the 127 TCGA genes onto the largest component - PPI_largest_comp
#         We get N genes out of 127 in the largest component  

driver_genes = open("gene_list.txt", "r")
tcga_genes = []
for gene in driver_genes:
	gene = gene.rstrip()
	tcga_genes.append(gene)
	
mapped_genes = []
for gene in tcga_genes:
	if gene in PPI_largest_comp.nodes():
		mapped_genes.append(gene)

	
print "The number of TCGA genes which map on the largest component of the PPI is", len(mapped_genes)	

#-------------------------------------------------------------------------------------------------------------#
# Step-4: Calculate average shortest path of these N genes in the largest component

def bfs_sp(gr, s): 
	""" Finds the shortest number of hops required to reach a node from s.
	Returns a dictionary with mapping: destination node from s -> no. of hops
	"""

	dist = {}
	q = deque([s])
	nodes_explored = set([s])
        for n in gr.nodes():
            if n == s: dist[n] = 0
        while len(q) != 0:
            node = q.popleft()
            for each in gr.neighbors(node):
                if each not in nodes_explored:
                    nodes_explored.add(each)
                    q.append(each)
                    dist[each] = dist[node] + 1
        return dist



list_distance1= []	
for n in mapped_genes:
	distance1 = bfs_sp(PPI_largest_comp, n)
	for k,v in distance1.items():
		if k in mapped_genes:
			list_distance1.append(v)
		
	

L1 = len(list_distance1)
print L1
y = 0
while True:
	try:
		list_distance1.remove(0)
	except:
		break	

L2  = len(list_distance1)
print L2
S = sum(list_distance1)
ASD = float (S)/L2
print "The average shortest distance of the 85 tcga genes mapped on to the largest component of the Reactome Interaction graph is", round(ASD, 2)

#------------------------------------------------------------------------------------------------------------#

# Step -5: Randomly pick up N = 85 genes from the largest component, 
#          Calculate the average shortest path among these randomly picked-up genes. 
#          Repeat this step 100 times to get a list of 100 values of average shortest paths.

nodes_lc = nx.nodes(PPI_largest_comp)

p = []
i = 0
while i < 5:
	random_sample2 = random.sample(nodes_lc, 85)
	list_distance2= []	
	for n in random_sample2:
		distance2 = bfs_sp(PPI_largest_comp, n)
		for k,v in distance2.items(): 
			if k in random_sample2:
				list_distance2.append(v)
	y = 0
	while True:
		try:
			list_distance2.remove(0)
		except:
			break			
	
	
	L4  = len(list_distance2)
	S2 = sum(list_distance2)
	asd = float (S2)/L4 
	p.append(round(asd, 2))
	i += 1

print p

stop = timeit.default_timer()

print "Time taken to run the code so far", (stop - start) 

#------------------------------------------------------------------------------------------------------------#

# Step-6: Use 100 values from Step-5 as the distribution to calculate p-value for the value from Step-4

# I ran the code with 100 iterations of Step-5 to get all the values in list_asd for this step
# It took 425.75 seconds 
# Probably not very efficient coding!


list_asd = [4.21, 4.28, 4.03, 4.16, 4.14, 4.24, 4.36, 4.05, 4.04, 4.1, 
			4.11, 4.27, 4.27, 4.11, 4.06, 4.1, 4.09, 3.94, 4.15, 3.95, 
			4.39, 3.97, 4.13, 4.26, 4.25, 4.3, 4.16, 4.07, 4.09, 4.15, 
			4.24, 4.21, 4.11, 4.14, 4.18, 4.3, 4.23, 4.06, 4.04, 4.15, 
			4.15, 4.39, 4.16, 4.14, 4.11, 4.13, 4.24, 4.16, 4.4, 4.14, 
			4.14, 4.15, 4.05, 4.05, 4.29, 4.16, 4.14, 4.03, 4.05, 4.17, 
			4.33, 4.18, 4.13, 4.11, 4.14, 4.14, 4.22, 4.35, 4.59, 4.17, 
			4.38, 4.13, 4.16, 4.07, 4.2, 4.2, 4.23, 4.28, 4.2, 4.18, 
			4.11, 4.0, 4.15, 4.09, 4.07, 4.15, 4.17, 4.23, 4.15, 4.19, 
			4.3, 4.29, 4.07, 4.17, 4.03, 4.22, 4.25, 4.26, 4.2, 4.04]

len_list_asd = len(list_asd)
sum_list_asd = sum(list_asd)

avg_list_asd = sum_list_asd/len_list_asd
# print avg_list_asd


sem = statistics.stdev(list_asd)/math.sqrt(len_list_asd)

# print sem

ucl = avg_list_asd + 1.98 * sem
lcl = avg_list_asd - 1.98 * sem

print "The 95-percent CI for the mean shortest distance of a graph which has randomly picked 85 nodes is", round(lcl,2),"-", round(ucl,2)
print "Since the value of the average shortest distance of the mapped 85 tcga genes does not lie in the above interval it is significantly different than the average shortest distance of a random graph."

