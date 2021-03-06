# SLEDGE: Synthetic LinEage Dataset GEnerator

Sledge is a python module which generates gene expression profiles mimicking scRNA-seq datasets.

The profiles are calculated according to the following principles:

* The genes are grouped into independent sets of modules.
* A module contains genes upregulated in a specific cell state.
* Each set of gene modules is associated with a specific type of trajectories.
* The trajectories can be either:
    * linear: gene expression is cascading in-between a linear sequence of transitions between cell states (e.g. A->B, B->C between cell states A, B and C).
    * bifurcating: a state can be followed by two child states (e.g. A->B, A->C).
    * cycling: the genes are expressed along a cycle of states (A->B, B->C, C->A)
* The depth in the trajectory associated with the first set of modules determines the timepoint of the cells.
* Each gene profile follow a normal distribution along its trajectory.
* A gaussian noise term can be added to gene expression.
* An extra set of unexpressed genes can also be joined.

<p align="center"><img src="schematic/schematic.png" width="100%"></p>
<p align="center"><i>The synthetic dataset is build upon 4 sets of gene modules: the first set is associated with the bifurcating cell lineage (21 states/modules), the second set is related to a simpler process with two linear transitions (3 states/modules), the third is mapped to a cycling set of 4 states (e.g cell cycle phases), and the fourth to a group of unexpressed genes. The asterisk relates how gene expression is represented on the heatmap (module 3, left) and along the differentation tree (level indicated in green, right).</i></p>

The synthetic dataset is stored as a set of csv files written in the working directory. The output folder also contains graph plots describing the cell state trajectories and heatmaps representing the gene expression matrix.

## How to Install

Sledge can be installed with `pip`.

```
pip install --upgrade https://github.com/juliendelile/Sledge/tarball/master
```

## How to Run

The dataset is generated with a single function.

```python
from sledge import generate_dataset

generate_dataset(
    n_cells                  = 500,
    n_tree_states            = [21, 3],
    p_branching              = [0.3, 0],
    n_genes_per_tree_state   = [200, 500],
    n_cc_states              = 4,
    n_genes_per_cc_state     = 500,
    n_unexpressed_genes      = 1000)
```

`generate_dataset` takes the following arguments:

| Argument  | Type | Description |
| --- | --- | --- |
| n_cells | int | the number of cells |
| n_tree_states | list of int | each integer determines the number of states composing a specific non-cycling trajectory. | 
| p_branching | list of float | the probability that a (non-cyclic) state is bifurcating. The sequence of transition is linear if equal to 0. |
| n_genes_per_tree_state | list of int | the number of genes associated with each state of a non-cycling trajectory. |
| n_cc_states | int | the number of cell cycle states. |
| n_genes_per_cc_state | int | the number of genes associated with each state of the cell cycle. |
| n_unexpressed_genes | int | the number of unexpressed genes. |
| noise_intensity | float | weight scaling the amplitude of the gaussian noise added to the gene expression levels. |
| seed | int | seed of the random number generator. |
| plot_size_factor | float | scaling factor for the nodes of the graph plots. |
