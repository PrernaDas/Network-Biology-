# Network-Biology-
The code uses the node and edge list from a simplified Reactome Interaction Network to build a Reactome Interaction Network Graph which has many components. It then takes the largest component and maps the 127 TCGA genes onto the largest component. It calculates the average shortest path between the genes that map on the largest component. It then randomly picks up the genes equivalent to those out of 127 TCGA that match and finds the average shortest path among these randomly picked-up genes. This step of randomly picking up genes and finding the average shortest path is repeated 100 times to get a list of 100 shortest paths. The key idea is that the average shortest path between the genes which are part of a network/pathway will be shorter than an equivalent number of random genes.
