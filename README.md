# PPRviz Reproducibility

## Single-level graph visualization

### Requirements
* Ubuntu
* python2.7 and python3.7
* sklearn, fa2l, networkx, tulip and networkit
### Location
* source code: `PPRviz-reproducibility/PPRVizS/`
* visualization output: `PPRviz-reproducibility/pprvizs_output/`
### Run layout algorithm
```
python benchmark.py [--data file_no] [--repeat times] [--mode support_mode] [--algo layout];
```
For example,
```
python2.7 benchmark.py --data 0 --repeat 10 --mode metrics --algo pprvizs;
```
### Input arguments
#### Datasets
Location: `PPRviz-reproducibility/dataset`

```
--data={0,1,2,3,4,5}
where
{0: "twego", 1: "fbego", 2: "wikiedit", 3: "physic", 4: "trust", 5: "scinet"}
```
#### Support algorithm
```
--algo={"pprvizs","ll","fa2","fr","mds","pivot_comp","simrank_comp","maxent","gf","le","lle","node2vec","sdne"}
```
##### Run in python2.7:
* `pprvizs`: our proposal PPRviz.
* `ll` and `fa2`: LinLog and ForceAtlas are supported by library `fa2l`
* `fr`: Fruchterman-Reingold is supported by library `networkx`
* `mds`: classical-MDS is supported by library `sklearn`
##### Run in python3.7:
* `pivot_comp`: pivot-MDS is supported by `tulip`
* `simrank_comp`: Simrank is supported by `networkx`
* `maxent`: Kadraw's single level is supported by `networkit`
* `gf,le,lle,node2vec,sdne`: node embedding methods are supported by [GEM](https://github.com/palash1992/GEM). We store the embeddings in `PPRviz-reproducibility/gem_pos/`.
#### Support mode
```
--mode={"plot","metrics"}
```
### benchmarking demo scripts
Visualization of all methods
```
bash run_plot.sh 
```
Metrics of all methods
```
bash run_metrics.sh 
```
## Multi-level graph visualization

### Requirement
* Ubuntu
* C++ 11
* Boost
* OpenMP
* google sparsehash
* Eigen
* cmake
### Dataset
`{6: "amazon", 7: "youtube", 8: "dblp", 9: "orkut", 10: "it", 11: "tw"}`
                             
### Clustering algorithm: Louvain+
#### Location
* source code: `PPRviz-reproducibility/louvain/`
* supergraph hierarchy output: `PPRviz-reproducibility/hierachy-output/` 
* leafnode partition output: `PPRviz-reproducibility/mapping-output/` 
#### Data format
Convert edgelist format dataset in `PPRviz-reproducibility/dataset/` into binary format for Louvain+ algorithm:
```
./convert -i ../dataset/amazon.txt -o ../dataset/amazon.bin;
```
#### Compile
```
cmake -DCMAKE_BUILD_TYPE=Release .; make;
```
#### Run clustering algorithm
```
./louvainplus [-f file_no] [-a algorithm] [-k partition_size] [-v] [-o];
```
where
```
-f: input file.
-a: 0 is conventional Louvain; 1 is Louvain+.
-k: threshold of partition size; default=25.
-v: verbose mode.
-o: output the partition.
```
For example
```
./louvainplus -f 6 -a 1 -k 25 -v;
```
### DPPR approximation algorithm: Tau-Push
#### Location
* source code: `PPRviz-reproducibility/approx_dnppr/`
* position matrix output: `PPRviz-reproducibility/[filename]_idx/`
#### Compile
```
cmake -DCMAKE_BUILD_TYPE=Release .; make;
```
#### Run DPPR approximation
```
./approx_dnppr [-f file_no] [-alg algorithm] [-build build_mode] [-sample query_times] [-random query_mode] [-embed compute_embedding];
```
where
```
-f: input file.
-alg: powiter; fora; foratp; fpsn; taupush.
-build: 1 is the mode of index construction and 0 otherwise.
-sample: number of queries.
-random: 1 is the random query and 0 is the largest-DPPR query.
-embed: 1 invokes stress majorization and 0 computes DPPR only.
```
#### Benchmarking demo scripts
Index construction of all methods
```
bash build_run.sh 
```
Random query of all methods
```
bash query_run.sh 
```
### Multi-level visualization demo
#### Location
* source code: `PPRviz-reproducibility/PPRVizl/`
* visualization output: `PPRviz-reproducibility/pprvizl_output/`
#### Visualize a multi-level path of amazon
```
python2.7 load-superppr-viz.py;
```
#### Result
​             ![Level-1](pprvizl_output/amazon-c0_l2_57.pdf) ![Level-0](pprvizl_output/amazon-c0_l1_3715.pdf)
