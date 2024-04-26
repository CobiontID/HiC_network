
# Extracting genomes from mixed samples with chromatin network embeddings

## Tools
### [Identifying groups of connected scaffolds](./node2vec)
Using information about pairwise chromatin connections between assembled sequences, this tool constructs a graph and allows interconnected groups of sequences to be clustered, labelled and visualised. Distinct genomes from contaminated samples can therefore be separated, even when reference-based taxonomic assignment is difficult. The procedure is based on the [node2vec](https://arxiv.org/abs/1607.00653) algorithm.

### [Hi-C map viewer](./map_viewer/)
Selectively display chromatin contact maps based on a list of contigs or scaffolds of interest, using a cooler file as input. Useful for displaying small fragments, or clusters identified by node2vec.


