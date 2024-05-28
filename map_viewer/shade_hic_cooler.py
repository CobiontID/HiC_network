#!/usr/bin/env python
# cw21@sanger.ac.uk

import argparse

import numpy as np
from scipy import sparse

import cooler
from iced import normalization

from holoviews.operation.datashader import datashade, rasterize, dynspread
import holoviews as hv
import matplotlib.pyplot as plt

hv.extension('bokeh')
import colorcet as cc
from holoviews import streams
from bokeh.resources import INLINE


#%%

parser = argparse.ArgumentParser(
    description='Generate a datashaded Hi-C map from a specified list of contigs or scaffolds')
parser.add_argument("--cooler", help="Cooler file to extract data from", required=True)
parser.add_argument("--contigs", help="Names or indices of contigs to stack, in order in which they are to be stacked.", required=False)
parser.add_argument("--contig_names", help="Indicate whether contig list contains names (otherwise, assume indices). Default True", type=bool, default=True)
parser.add_argument("--sort", help="Sort contigs by number", type=bool, default=True)
parser.add_argument("--name", help="Name for output", required=True)
parser.add_argument("--palette", help="Pick colour palette", default="blues")


args = parser.parse_args()

print(args)

#required
cool = args.cooler
contigfile = args.contigs
contig_use_names = args.contig_names
name = args.name
palette = args.palette

sort_contigs = args.sort

#%%
def load_from_cool(clr, s1, s2):
    """Load data for scaffolds s1 and s2"""
    a = clr.matrix(balance=False)[chromstarts[s1]:chromstarts[s1+1],chromstarts[s2]:chromstarts[s2+1]]
    return a

def norm_matrix(sparse_matrix):
    """Convert sparse matrix to dense and normalise"""
    a_norm = np.asarray(normalization.ICE_normalization(sparse_matrix).todense())
    a_norm[np.isnan(a_norm)] = 0.
    return a_norm

def stack_contigs(clr, contigs):
    """Extract sections of the cooler matrix, stack, and normalise.
    Return matrix and lengths of each contig's corresponding segment"""
    s = {}
    for i in range(len(contigs[:])):
        s[i] = np.hstack([load_from_cool(clr, contigs[i], j) for j in contigs])
    
    m = np.vstack([s[i] for i in range(len(contigs))])
    sparse_matrix = sparse.coo_matrix(m)
    lengths = [s[i].shape[0] for i in range(len(contigs))]
    return norm_matrix(sparse_matrix), lengths


# Sum connections
# squashed =  np.array([[np.sum(load_from_cool(clr, contigs[i], j)) for j in contigs] for i in range(len(contigs))])

# Load cooler
clr = cooler.Cooler(cool)

# Get contig ids
if contig_use_names is True:
    contig_names = np.loadtxt(contigfile, dtype="str")
    contigs = np.array([clr.chromnames.index(i) for i in contig_names])
else:
    contigs = np.loadtxt(contigfile)

if sort_contigs is True:
    contigs = np.sort(contigs)

chromstarts = [clr.extent(i)[0] for i in clr.chromnames]

# Load and stack data
a_norm, lengths = stack_contigs(clr, contigs)

# Plotting
ticks = [z for z in zip([sum(lengths[:i+1]) for i,l in enumerate(lengths)],contigs+1)][:]
#grid = {'x': [0, ticks[0][0]], 'y': [ticks[0][0], 0], 'x1': [len(a_norm), ticks[0][0]], 'y1': [ticks[0][0], len(a_norm)]}

grid = {'x': [i[0] for i in ticks] + [0]*len(ticks),
        'y': [0]*len(ticks) + [i[0] for i in ticks],
        'x1': [i[0] for i in ticks] + [len(a_norm)]*len(ticks),
        'y1': [len(a_norm)]*len(ticks) + [i[0] for i in ticks]
        }

shaded = datashade(hv.Image(np.flipud(a_norm), bounds=(0,0,len(a_norm), len(a_norm))), cmap=palette).opts(
    width=1000, height=1000, invert_yaxis=True, xticks=ticks, yticks=ticks, xlabel="", ylabel="") * hv.Segments(
        grid, kdims=['x', 'y', 'x1', 'y1']).opts(color="grey", line_width=0.5)

print(f"Saving plot as {name}.html")
hv.save(shaded, f'{name}.html')
