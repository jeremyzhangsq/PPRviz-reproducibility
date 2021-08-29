import numpy as np
import networkx as nx
from PPRVizS.pprviz import rotate


def d(u, v):
    return np.sqrt(np.sum((u - v) ** 2, axis=0))


"""
The sum of node replusion closys in the normalized drawing over all pairs of nodes.
The node replusion from ant pair of nodes is the reciprocal of the square of the distance between them.
The lower the better
"""

# See Section 2.2.1 in https://crpit.scem.westernsydney.edu.au/confpapers/CRPITV60Lee.pdf
# Or Section 3 in https://dl.acm.org/doi/pdf/10.1145/234535.234538
def get_node_distribution(normX):
    n = len(normX)
    cost = 0
    for u in range(n):
        for v in range(n):
            if u <= v:
                continue
            dist = d(normX[u], normX[v])
            if dist==0:
                cost = float('Inf')
                break
            else:
                cost += 1 / (dist ** 2)
    return cost




# see Edge length variation Section 4.1.4 in https://arxiv.org/pdf/1710.04328.pdf
def uniform_edge_length_coefficient_variance(edges, X):
    length = []
    for e in edges:
        length.append(d(X[e[0]], X[e[1]]))

    lsigma = np.std(length)
    lmu = np.mean(length)
    lcv = lsigma / lmu
    return lcv


def best_fit_layout(X):
    rotatedX = rotate(
        X)  # rotate to fix minimal bounding box, i.e., the rectangular covering all points with minimal area
    width = rotatedX[:, 0].max() - rotatedX[:, 0].min()
    height = rotatedX[:, 1].max() - rotatedX[:, 1].min()

    wscale = width / 4  # fix the width of the canvas to be 4
    hscale = width / 4  # scale height as width does

    normX = rotatedX
    normX[:, 0] = normX[:, 0] / wscale
    normX[:, 1] = normX[:, 1] / hscale
    return normX



def eva(G, X):
    edges = G.edges()
    try: normX = best_fit_layout(X)
    except: normX = X
    try: a = get_node_distribution(normX)
    except: a = np.infty
    b = uniform_edge_length_coefficient_variance(edges, normX)
    return a,b


if __name__ == '__main__':
    G = nx.karate_club_graph()
    n, m = G.number_of_nodes(), G.number_of_edges()
    print(n, m)
    edges = G.edges()
    X = np.random.uniform(low=0, high=np.sqrt(2 * n - 2), size=(n, 2))
    eva(G, X)
