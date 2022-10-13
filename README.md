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
python benchmark.py [--data file_name] [--repeat times] [--mode support_mode] [--alg layout];
```
For example,
```
python2.7 benchmark.py --data twego --repeat 10 --mode metrics --alg pprvizs;
```
### Input arguments

#### Datasets

Location: `PPRviz-reproducibility/dataset`
```
--data={"twego","fbego","wikiedit","physic","trust","scinet"}
```

#### Support algorithm

```
--alg={"pprvizs","ll","fa2","fr","mds","pivot_comp","simrank_comp","maxent","gf","le","lle","node2vec","sdne"}
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

`{"amazon", "youtube", "dblp", "orkut", "it", "tw"}`                      

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
./louvainplus [-f file_name] [-alg algorithm] [-k partition_size] [-s random seed] [-v 0/1] [-o 0/1];
```
where
```
-f: input file.
-a: 0 is conventional Louvain; 1 is the adapted version called Louvain+.
-k: threshold of partition size; default=25.
-s: random seed.
-v: 1 is verbose mode.
-o: 1 is for output the partition.
```
For example
```
./louvainplus -f amazon -a 1 -k 25;
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
./approx_dnppr [-f file_no] [-alg algorithm] [-k size_limit] [-build build_mode] [-sample query_times] [-random query_mode] [-embed compute_embedding];
```
where
```
-f: input file.
-alg: powiter; fora; foratp; fpsn; taupush.
-k: cluster size limit (default=25)
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

#### Evaluate the multi-level result

* compute the metric scores for PPRviz variants
```
python2.7 load-superppr-viz.py --data amazon --k 25 --alg taupush --mode metrics;
```


* visualize the high-level results for PPRviz variants and the outputs are stored as [here](./pprvizl_output).
```
python2.7 load-superppr-viz.py --data amazon --k 25 --alg taupush --mode plot;
```

