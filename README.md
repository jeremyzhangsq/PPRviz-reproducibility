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

### more metrics evaluation for author feedback

We show the CR and AR scores in CR/AR format.
Note that a lower score indicates a better quality.

```
+-----------+-----------------+------------------+-----------------+------------------+------------------+------------------+------------------+------------------+-------------------+-------------------+------------------+------------------+
|           |      PPRviz     |        FR        |    ForceAtlas   |      LinLog      |       CMDS       |       PMDS       |      SimRank     |      GFactor     |        SDNE       |       LapEig      |        LLE       |     Node2vec     |
+-----------+-----------------+------------------+-----------------+------------------+------------------+------------------+------------------+------------------+-------------------+-------------------+------------------+------------------+
|   TwEgo   |  1.00E+00/0.00  |   5.10E+00/2.60  |  9.00E-01/0.26  |   9.00E-01/0.65  |   4.00E+00/0.00  |   2.00E+00/5.00  |   4.00E+01/3.45  |  5.70E+01/13.97  |   2.24E+02/59.99  |   7.50E+01/25.08  |   6.60E+01/3.87  |  7.80E+01/12.72  |
+-----------+-----------------+------------------+-----------------+------------------+------------------+------------------+------------------+------------------+-------------------+-------------------+------------------+------------------+
|   FbEgo   |  3.10E+02/40.89 |  3.75E+02/58.44  |  2.59E+02/52.27 |  3.25E+02/63.70  |  2.75E+02/45.41  |  3.02E+02/98.12  |  4.37E+02/64.34  |  1.42E+03/271.35 |  7.04E+03/958.00  |  2.09E+03/643.88  |  1.86E+03/741.00 |  2.16E+03/296.71 |
+-----------+-----------------+------------------+-----------------+------------------+------------------+------------------+------------------+------------------+-------------------+-------------------+------------------+------------------+
|  Wiki-ii  | 1.12E+03/408.36 |  1.22E+03/402.04 | 6.31E+02/521.96 |  6.92E+02/478.43 |  1.47E+03/514.18 | 2.46E+03/3147.44 | 4.86E+03/1457.52 |  8.59E+03/614.29 |  3.09E+04/4495.00 |  1.12E+04/2695.72 | 9.33E+03/2251.41 | 1.03E+04/1012.85 |
+-----------+-----------------+------------------+-----------------+------------------+------------------+------------------+------------------+------------------+-------------------+-------------------+------------------+------------------+
| Physician | 9.74E+03/459.70 |  8.80E+03/894.17 | 6.65E+03/391.67 |  6.39E+03/440.03 |  1.86E+04/419.35 |  8.18E+03/582.90 |  3.05E+04/781.61 | 4.63E+04/1376.33 |  1.83E+05/5588.83 |  3.07E+04/4739.49 | 9.45E+04/5339.00 | 2.70E+04/2404.98 |
+-----------+-----------------+------------------+-----------------+------------------+------------------+------------------+------------------+------------------+-------------------+-------------------+------------------+------------------+
| FilmTrust | 1.12E+04/671.64 |  1.30E+04/827.74 | 6.32E+03/782.83 |  5.54E+03/758.24 | 1.95E+04/1147.72 | 1.26E+04/2016.99 | 1.52E+04/1018.14 | 1.56E+05/1510.53 |  1.40E+05/7794.39 |  1.32E+05/6627.20 | 1.47E+05/3725.22 | 1.21E+05/3438.33 |
+-----------+-----------------+------------------+-----------------+------------------+------------------+------------------+------------------+------------------+-------------------+-------------------+------------------+------------------+
|   SciNet  | 1.20E+04/751.06 | 6.73E+03/1421.84 | 5.23E+031301.16 | 5.22E+03/1517.37 | 3.41E+05/1204.97 | 1.50E+04/7520.11 | 1.87E+04/1097.92 | 4.17E+05/3603.67 | 3.66E+06/16242.00 | 2.38E+05/10037.90 |        -/-       | 2.34E+05/5167.96 |
+-----------+-----------------+------------------+-----------------+------------------+------------------+------------------+------------------+------------------+-------------------+-------------------+------------------+------------------+
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

#### Evaluate the multi-level result

* compute the metric scores for PPRviz variants
```
python2.7 load-superppr-viz.py --data 5 --mode metrics;
```

The ND/ULCv/CR/AR scores on the top level of SciNet is listed as follows.
```
Exact
8.77E+01/0.87/1.20E+01/67.65
FORA-TP
8.96E+01/0.88/9.00E+00/67.57
Tau-Push
4.57E+01/0.86/1.00E+01/66.27
```

* visualize the high-level results for PPRviz variants and the outputs are stored as [here](./pprvizl_output).
```
python2.7 load-superppr-viz.py --data 5 --mode plot;
```

* the top-level results in SciNet for PPRviz with Exact, FORA-TP and Tau-Push are shown as follows:

![Exact](./pprvizl_output/scinet-c0_l2_0-powiter.pdf)

![FORA-TP](./pprvizl_output/scinet-c0_l2_0-foratp.pdf)

![Tau-Push](./pprvizl_output/scinet-c0_l2_0-taupush.pdf)
