import sys
sys.path.insert(0, '../')
sys.path.insert(0, './')
from fa2l import force_atlas2_layout
import matplotlib
matplotlib.use('Agg')
from PPRVizS.evaluate import eva
from PPRVizS.pprviz import *
# from PPRVizS.config import filelist
import argparse
if sys.version_info >= (3, 0):
    from tulip import tlp
    import networkit as nt


# pip install tulip-python
# https://tulip.labri.fr/Documentation/current/tulip-python/html/tulippluginsdocumentation.html#algorithmpluginsdoc
def save_pivot_mds(G):
    nodesDict = [None] * len(G)
    graph = tlp.newGraph()
    for node in G.nodes():
        tnode = graph.addNode()
        nodesDict[node] = tnode

    for (u, v) in G.edges():
        graph.addEdge(nodesDict[u], nodesDict[v])


    start = time.time()
    params = tlp.getDefaultPluginParameters('Pivot MDS (OGDF)', graph)

    # either create or get a layout property from the graph to store the result of the algorithm
    resultLayout = graph.getLayoutProperty('resultLayout')
    successFlag = graph.applyLayoutAlgorithm('Pivot MDS (OGDF)', resultLayout, params)
    # or store the result of the algorithm in the default Tulip layout property named 'viewLayout'
    # success = graph.applyLayoutAlgorithm('Pivot MDS (OGDF)', params)
    t = time.time() - start

    pos = np.zeros((len(G), 2))
    for i in range(len(nodesDict)):
        tnode = nodesDict[i]
        (x, y, z) = resultLayout.getNodeValue(tnode)
        pos[i, 0] = x
        pos[i, 1] = y
    f = "../gem_pos/{}_{}.txt".format(args.data, "pivot")
    np.save(f, pos)
    exit(0)

def lin_log(G):
    positions = force_atlas2_layout(G, iterations=700, gravity=0,
                                    barnes_hut_theta=1.2, barnes_hut_optimize=True,
                                    edge_weight_influence=1, jitter_tolerance=1,
                                    scaling_ratio=2, strong_gravity_mode=False,
                                    lin_log_mode=True)
    return np.array([positions[i] for i in range(n)])


def force_atlas(G):
    positions = force_atlas2_layout(G, iterations=700, gravity=0,
                                    barnes_hut_theta=1.2, barnes_hut_optimize=True,
                                    edge_weight_influence=1, jitter_tolerance=1,
                                    scaling_ratio=2, strong_gravity_mode=False,
                                    lin_log_mode=False)

    return np.array([positions[i] for i in range(n)])


def classical_mds(G):
    _,pos = dist_mds(G=G)
    return pos

def pprviz_s(G):
    return pprdegree_relativemds(G=G,is_vanilla=False)


def save_simrankviz_vanilla(G):
    sim = nx.simrank_similarity(G)
    lol = [[sim[u][v] for v in sorted(sim[u])] for u in sorted(sim)]
    simRank = csr_matrix(lol)
    X = pprdegree_relativemds(PPR=simRank, is_vanilla=True, is_vanilla_scale=False)
    f = "../gem_pos/{}_{}.txt".format(args.data, "simrank")
    np.save(f, X)
    exit(0)


def maxent_stress(G):
    Gnew = nt.Graph(n)
    for u,v in G.edges():
        Gnew.addEdge(u,v)
    # nt.overview(Gnew)
    init = save_pivot_mds(G)
    # init = np.random.uniform(size=(n,2))
    maxent = nt.viz.MaxentStress(Gnew,2,1,init,0.001)
    maxent.run()
    maxent.scaleLayout()
    pos = maxent.getCoordinates()
    return np.array(pos)

def fruchterman_reingold(G):
    out = nx.spring_layout(G=G, dim=2)
    pos = np.zeros((n,2))
    for i in out:
        pos[i, 0] = out[i][0]
        pos[i, 1] = out[i][1]
    return pos

def run(alg):
    if alg == "pprvizs":
        X = pprviz_s(G)
    elif alg == "simrank_comp":
        X = save_simrankviz_vanilla(G)
    elif alg == "pivot_comp":
        X = save_pivot_mds(G)
    elif alg == "maxent":
        X = maxent_stress(G)
    elif alg == "mds":
        X = classical_mds(G)
    elif alg == "fa2":
        X = force_atlas(G)
    elif alg == "ll":
        X = lin_log(G)
    elif alg == "fr":
        X = fruchterman_reingold(G)
    elif alg in ["gf","le","lle","node2vec","sdne","simrank","pivot"]:
        X = np.load("../gem_pos/{}_{}.txt.npy".format(args.data,alg))
    else:
        print("error input algorithm:(")
        exit(-1)
    return X

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process...')
    parser.add_argument('--data', type=str, default="twego", help='graph dataset id')
    parser.add_argument('--repeat', type=int, default=1, help='repeating number')
    parser.add_argument('--algo', type=str, default="pprvizs", help='repeating number')
    parser.add_argument('--mode',type=str,default="metrics",help='plot or metrics')
    args = parser.parse_args()
    path = "../dataset/" + args.data + ".txt"
    G = nx.read_edgelist(path, nodetype=int)
    n, m = G.number_of_nodes(), G.number_of_edges()
    edges = G.edges()
    if args.mode == "metrics":
        num_repeat = args.repeat
        nd = np.zeros(num_repeat)
        ulcv = np.zeros(num_repeat)
        ar = np.zeros(num_repeat)
        try:
            for i in range(num_repeat):
                start = time.time()
                X = run(args.algo)
                nd[i],ulcv[i],ar[i] = eva(G, X)
        except Exception as e:
            print(e)
            print("NaN NaN NaN")
        else:
            print ("{:.2E}/{:.2f}/{:.2f}".format(nd.mean(),ulcv.mean(),ar.mean()))
            # output visualization layout
    if args.mode == "plot":
        print("plotting {} on {}".format(args.algo, args.data))
        try:
            X = run(args.algo)
        except Exception as e:
            print(e)
        else:
            is_scale = True
            if args.algo == "maxent":
                is_scale = False
            plot(pos=X, edges=edges, scale=is_scale,
                     name="../pprvizs_output/{}-{}".format(args.data, args.algo))
