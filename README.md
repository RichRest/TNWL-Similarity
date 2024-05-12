# TNWL

This project contains an implementation of the Temporal Neighborhood Weisfeiler-Lehman algorithm for temporal graphs. Included are a temporal graph kernel, a similarity measure based on that and a label analysis tool for two differing graphs.
Also, with the help of the Voro++ library a transformation of chemical simulations to their temporal Delaunay graph is implemented.

## Requirements

- cmake 3.10
- Basic libraries
- Voro++

## Input

The folder containing any dataset should be named as `dataset` and be stored in the folder "datasets". The source code supports the following three data formats:

### XYZ
This input is used for chemical simulations. The file must be stored as `dataset.xyz`

Every time step follows the format
```
   n
*line of chemical specifications*
T x y z ...
T x y z ...
...
```
with `n` denoting the number of atoms and `T` denoting the atom type.

### TU
Fits the data from the TUDataset https://chrsmrrs.github.io/datasets/. The files `dataset_A.txt`, `dataset_edge_attributes.txt`, `dataset_graph_indicator.txt`, and `dataset_node_labels.txt` are necessary. Note that the graphs are stored in a single graph and the partition based on `dataset_graph_indicator.txt` is realized as edge labels.

### vwt5
This input can be used for datasets with very large vertex labels as these will be compressed. The file must be stored as `dataset.txt`.

The file is a list of all edges in the format
```
v w t ... l
```
with `v,w` denoting the edge and tail, `t` the time label, and `l` the label of the edge. The label has to be in the fifth column.

## Main functionalities

`create_TG(...)`: Creates a temporal graph from the given dataset with the specified input format

`TNWL(...)`: Performs the TNWL algorithm (as stated in Section 3.2) on a list of graphs

`label_analysis(...)`: Performs the label analysis stated in Section 4.1 and outputs plots in the folders `event_`

`print_gram_of_events(...)`: Takes a graph initialized from the TUDataset and computes the gram matrix of the pairwise TNWL similarities. The resulting gram matrix is already normalized.
