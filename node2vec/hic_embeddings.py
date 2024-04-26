# %%
import os.path

import numpy as np
import matplotlib.pyplot as plt
import networkx

import smart_open
smart_open.open = smart_open.smart_open
from gensim.models import Word2Vec
import random
import node2vec


#import datashader as ds
#from holoviews.operation import decimate
from holoviews.operation.datashader import datashade, rasterize, dynspread
from holoviews import dim, opts
import holoviews as hv
import matplotlib.pyplot as plt
import pandas as pd

#from sklearn.preprocessing import MinMaxScaler
import numpy as np
hv.extension('bokeh')
#import colorcet as cc
import hvplot.networkx as hvnx

from sklearn.cluster import KMeans

from scipy.linalg import fractional_matrix_power


import argparse

##
#%%
parser = argparse.ArgumentParser(
    description='Plot decomposed contigs')
parser.add_argument("--infile", help="input .npy file of connections", required="True")
parser.add_argument("--outfile", help="Target file for plot, must end in .html", required=True) #FIXME
parser.add_argument("--name", help="Name for outputs", required=True)
parser.add_argument("--k", help="Number of clusters", default=2, type=int)
parser.add_argument("--dimensions", help="Dimensions for embeddings", default=128, required=False, type=int)
parser.add_argument("--save_embeddings", help="Whether to save embeddings", default="F", required=False)
parser.add_argument("--p", help="p parameter", default=1, type=float, required=False)
parser.add_argument("--q", help="q parameter", default=1, type=float, required=False)


args = parser.parse_args()
print(args)

infile = args.infile
outfile = args.outfile
n_clusters = args.k
name = args.name
save_embeddings = args.save_embeddings
dimensions = args.dimensions
p = args.p
q = args.q

# %%

remove_disc = "False"


# Load connectivity matrix, zero the diagonal (connections with self), and produce plot

a = np.load(infile)
# Check symmetry
try:
    assert np.all(a.T == a)
except AssertionError:
    print("Matrix is not symmetrical") 

# %%

np.fill_diagonal(a, 0)
plt.figure(figsize=(15,15))
plt.imshow(a > 0, cmap="gray")
plt.yticks(range(len(a)))
#plt.show()
plt.savefig("{}.matrix_yes_no.png".format(name))
plt.close()

# Remove unconnected scaffolds if specified
# Caution, this breaks numbering
if remove_disc == "True":
    ix = np.where(np.sum(a, axis=0) > 0)[0]
    print(ix)
    a = a[ix][:,ix]

G = networkx.from_numpy_matrix(a > 0)

# %%

# Run random walks with specified parameters
# Generate node embeddings with the specified dimensions using gensim's Word2Vec using the skip-gram algorithm

output = "run.walks"
dimensions = dimensions
window_size = 20
workers = 10
iter = 1

node2vec.run_model(G, output, dimensions, window_size, workers, iter, p, q)

# %%

X = np.loadtxt(output, skiprows=1) # load the embedding of the nodes of the graph
# sort the embedding based on node index in the first column in X
X=X[X[:,0].argsort()]
Z=X[0:X.shape[0],1:X.shape[1]] # remove the node index from X and save in Z

# If requested, save the sorted node embeddings as a numpy matrix
if save_embeddings == "T":
    np.save("{}.embeddings.npy".format(name), Z)

# Apply KMeans to embeddings, save cluster labels to file, and add node labels to network graph
kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(Z)
labels = kmeans.labels_ 

with open("{}.labels.txt".format(name), "w") as of:
    for l in labels:
        of.write(str(l)+"\n")

for i, _ in enumerate(G.nodes()):
    G.nodes[i]['type'] = labels[i]

#G.remove_nodes_from(list(networkx.isolates(G)))

# Plot network graph with labels
print("Rendering network graph.")
from holoviews.operation.datashader import datashade, bundle_graph

bundled = bundle_graph(hv.Graph.from_networkx(G, networkx.spring_layout))
node_labels = hv.Labels(bundled.nodes, ['x', 'y'], 'index')
hv.save((datashade(bundled) * bundled.nodes.opts(color="type", size=15)).opts(width=800,height=800), outfile)

#hv.save(hvnx.draw_spring(G, with_labels=True, edge_color="grey", node_color=labels, node_cmap=cc.glasbey_cool).opts(width=700, height=700), outfile) #"network_embeddings_fthualb1_recolour.html"

# Save datashaded matrix as html
shaded = datashade(hv.Image(a), cmap="viridis").opts(width=500, height=500)
hv.save(shaded, "{}.shaded_matrix.html".format(name))

# Re-sort the simplified connectivity matrix by assigned label
# Note: The order of the labels may not produce a coherent plot

idx = np.concatenate([np.where(np.array(labels) == k)[0] for k in range(n_clusters)])
plt.figure(figsize=(15,15))
plt.imshow(a[idx][:,idx] > 0, cmap="gray")
plt.yticks(range(len(a)))
#plt.show()
plt.savefig("{}.matrix_yes_no_sort.png".format(name))
plt.close()

# Plot 2D PCA of node2vec embeddings, colour nodes by assigned labels

print("Plotting 2D representation of embeddings.")

from sklearn.decomposition import PCA

pca = PCA(n_components=2)
transformed = pca.fit_transform(Z)

plt.figure(figsize=(10,10))
plt.scatter(transformed[:,0], transformed[:,1], c=labels)
plt.savefig("{}.embeddings_pca.png".format(name))
plt.close()
