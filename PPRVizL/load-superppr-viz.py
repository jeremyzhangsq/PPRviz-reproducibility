import matplotlib
import argparse
import sys
sys.path.insert(0, '../')
matplotlib.use('Agg')
from PPRVizS.pprviz import *
from PPRVizS.config import filelist
from PPRVizS.evaluate import eva, rescale_layout, best_fit_layout
from scipy import sparse

# global supernode_list, mapping_data, A, n
A=None
Gfull = None
D = None
supernode_list=dict()
mapping_data=dict()
PPRDeg = None
Xmds = None
r = None
super2id = dict()
id2super = []
component = 0
level = 0

def load_ppr(path, size):
    global PPRDeg
    with open(path+".src", "rb") as fin:
        row = np.fromfile(fin, dtype=np.int32)
    with open(path+".dst", "rb") as fin:
        col = np.fromfile(fin, dtype=np.int32)
    with open(path+".ppr", "rb") as fin:
        data = np.fromfile(fin, dtype=np.double)
    assert (len(row)==len(col))
    assert (len(row)==len(data))
    # print (len(row),len(col))
    row1 = [super2id[i] for i in row]
    col1 = [super2id[i] for i in col]
    PPRDeg = sparse.csr_matrix((data, (row1, col1)), shape=(size, size))

def load_position(path):
    global Xmds, r
    with open(path+".x", "rb") as fin:
        x = np.fromfile(fin, dtype=np.double)
    with open(path+".y", "rb") as fin:
        y = np.fromfile(fin, dtype=np.double)
    with open(path+".r", "rb") as fin:
        r = np.fromfile(fin, dtype=np.double)
    # print (len(row),len(col))
    Xmds = np.vstack([x,y])
    Xmds = Xmds.T

def load_community(hiefname, mapfname):
    global supernode_list, mapping_data
    with open(hiefname, "r") as fin:
        while True:
            l = fin.readline()
            if len(l) == 0:
                break
            name = l.rstrip("\r\n")
            l = fin.readline()
            size = int(l.rstrip("\r\n"))
            child = [int(fin.readline().rstrip("\r\n")) for i in range(size)]
            supernode_list[name] = child

    with open(mapfname, "r") as fin:
        while True:
            l = fin.readline()
            if len(l) == 0:
                break
            name = l.rstrip("\r\n")
            l = fin.readline()
            size = int(l.rstrip("\r\n"))
            child = [int(fin.readline().rstrip("\r\n")) for i in range(size)]
            mapping_data[name] = child


def split_supernode(name):
    global component, level
    lst = name.split("_")
    component = int(lst[0][1:])
    level = int(lst[1][1:])
    sid = int(lst[2])
    return sid

def get_children(node):
    global super2id,id2super
    split_supernode(node)
    if level == 1:
        childsize = len(mapping_data[node])
        id2super = mapping_data[node]
    else:
        childsize = len(supernode_list[node])
        id2super = supernode_list[node]
    super2id = {id2super[i]: i for i in range(childsize)}

    return childsize

def get_subgraph(cluster):
    s = len(cluster)
    cm1 = np.zeros((s, n))
    nodeweight = [1] * s

    for i in range(s):
        supernode = cluster[i]
        if supernode in mapping_data:
            cidx = mapping_data[supernode]
            nodeweight[i] = len(cidx)
        else:
            cidx = [supernode]
            nodeweight[i] = 1
        cm1[i, :] = np.sum(A[cidx, :], axis=0)

    As = np.zeros((s, s))
    for i in range(s):
        supernode = cluster[i]
        if supernode in mapping_data:
            cidx = mapping_data[supernode]
        else:
            cidx = [supernode]
        As[:, i] = np.sum(cm1[:, cidx], axis=1)

    subG = nx.from_scipy_sparse_matrix(sparse.csr_matrix(As))

    return subG, nodeweight


def viz(target,alg,k):
    start = time.time()
    if level==1:
        cluster = id2super
    else:
        cluster = ["c"+str(component)+"_l"+str(level-1)+"_"+str(i) for i in id2super]

    G, nodeweight = get_subgraph(cluster=cluster)
    edges = G.edges()
    t = time.time() - start
    print ("subgraph time:{}s n:{} m:{}".format(t,G.number_of_nodes(),G.number_of_edges()))

    start = time.time()

    t = time.time() - start
    print ("zoom-in time:{}s".format(t))

    plot(pos=Xmds, edges=edges, radius = r, scale=False,
         name="../pprvizl_output/{}-{}-{}-{}".format(filelist[dataid], target,alg,k))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process...')
    parser.add_argument('--data', type=int, default=5, help='graph dataset id')
    parser.add_argument('--mode', type=str, default="metrics", help='plot or metrics')
    args = parser.parse_args()
    dataid = args.data
    dataname = filelist[dataid]


    # print("loading clusters...")

    # print("loading edges...")
    # path = "/home/zhangsq/gviz-ppr/dataset/" + dataname +".txt"
    fpath = "../dataset/" + dataname + ".txt"
    Gfull = nx.read_edgelist(fpath, nodetype=int)
    A = nx.adjacency_matrix(Gfull)
    n = Gfull.number_of_nodes()

    algos = ["powiter","taupush"]


    for k in [15, 20, 25]:
        hiefname = '../louvain/hierachy-output/{}_{}.dat'.format(dataname,k)
        mapfname = '../louvain/mapping-output/{}_{}.dat'.format(dataname,k)
        load_community(hiefname, mapfname)
        storename = "../{}_idx/{}ds250_{}".format(dataname, dataname,k)
        zoompath = "c0_l2_0"
        childsize = get_children(zoompath)
        cluster = ["c0_l1_" + str(i) for i in id2super]

        G, nodeweight = get_subgraph(cluster=cluster)
        for algo in algos:
            print(algo)
            pospath = storename + "{}_{}".format(zoompath,algo)
            load_position(pospath)
            if args.mode == "plot":
                viz(zoompath,algo,k)
            elif args.mode == "metrics":
                nd,ulcv,cr,ar = eva(G, Xmds)
                print ("{:.2E}/{:.2f}/{:.2E}/{:.2f}".format(nd,ulcv,cr,ar))
            else:
                exit(-1)


