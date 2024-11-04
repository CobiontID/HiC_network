
# Extracting genomes from mixed samples with chromatin network embeddings

![image](https://github.com/CobiontID/HiC_network/assets/10507101/6ce87f4f-24f4-4849-b9d9-8951a152082a)


## Tools
### [Identifying groups of connected scaffolds](./node2vec)
Using information about pairwise chromatin connections between assembled sequences, this tool constructs a graph and allows interconnected groups of sequences to be clustered, labelled and visualised. Distinct genomes from contaminated samples can therefore be separated, even when reference-based taxonomic assignment is difficult. The procedure is based on the [node2vec](https://arxiv.org/abs/1607.00653) algorithm.

### [Hi-C map viewer](./map_viewer/)
Selectively display chromatin contact maps based on a list of contigs or scaffolds of interest, using a cooler file as input. Useful for displaying small fragments, or clusters identified by node2vec.

## Citation
If you use code from this repository, please cite: _Kudoa genomes from contaminated hosts reveal extensive gene order conservation and rapid sequence evolution_ https://www.biorxiv.org/content/10.1101/2024.11.01.621499v1

