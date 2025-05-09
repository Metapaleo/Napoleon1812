#!/usr/bin/env python3
import sys
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# Read input from stdin
data = [line.strip().split('\t') for line in sys.stdin]
header = data[0]
rows = data[1:]

# Extract positions and values
positions = [row[0] for row in rows]
matrix = np.array([list(map(int, row[1:])) for row in rows])

# Calculate pairwise distances (weighting 2s more than 1s)
def custom_dist(u, v):
    weights = np.where((u == 2) | (v == 2), 2, 1)
    return np.sum(weights * (u != v))

# Cluster rows using hierarchical clustering
pairwise_dists = pdist(matrix, metric=custom_dist)
linkage_matrix = linkage(pairwise_dists, method='average')
row_order = leaves_list(linkage_matrix)

# Print reordered table
print('\t'.join(header))
for i in row_order:
    print('\t'.join([positions[i]] + list(map(str, matrix[i]))))
