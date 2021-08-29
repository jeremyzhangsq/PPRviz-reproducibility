import numpy as np
import networkx as nx
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from sklearn import preprocessing
from scipy.sparse import identity, csr_matrix
from sklearn import manifold
import time
from PPRVizS.mds import relativeMDS
from PPRVizS.MinimumBoundingBox import MinimumBoundingBox
from math import cos, sin


def rotate(pos):
    points = []
    n = pos.shape[0]
    for u in range(n):
        points.append((pos[u, 0], pos[u, 1]))

    bounding_box = MinimumBoundingBox(points)  # returns namedtuple
    unit_vector_angle = bounding_box.unit_vector_angle
    newpos = pos
    for u in range(n):
        x = pos[u, 0] * cos(unit_vector_angle) + pos[u, 1] * sin(unit_vector_angle)
        y = pos[u, 1] * cos(unit_vector_angle) - pos[u, 0] * sin(unit_vector_angle)
        newpos[u, 0] = x
        newpos[u, 1] = y

    # for u in range(n):
    #    for v in range(n):
    #        print(d(pos[u],pos[v]), d(newpos[u],newpos[v]))
    return newpos


def get_boundary(pos, radius):
    xmin = np.Inf
    ymin = np.Inf
    xmax = np.NINF
    ymax = np.NINF
    for i in range(len(pos)):
        xminus = pos[i, 0] - radius[i]
        xplus = pos[i, 0] + radius[i]
        yminus = pos[i, 1] - radius[i]
        yplus = pos[i, 1] + radius[i]
        xmin = xminus if xminus < xmin else xmin
        xmax = xplus if xplus > xmax else xmax
        ymin = yminus if yminus < ymin else ymin
        ymax = yplus if yplus > ymax else ymax

    return xmin, xmax, ymin, ymax


def plot(pos, edges=None, radius=None, label=None, scale=True, name="ppr"):
    if scale == True:
        try:
            pos = rotate(pos)
        except:
            pass
    n = len(pos)

    if radius is not None:
        x_min, x_max, y_min, y_max = get_boundary(pos, radius)
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)

    plt.autoscale(False)
    fig, ax = plt.subplots()

    if radius is not None:
        height = plt.rcParams.get('figure.figsize')[1]
        y_range = pos[:, 1].max() - pos[:, 1].min()
        dpi = plt.rcParams.get('figure.dpi')
        scaler = 1.0 * dpi * height / y_range
        maxs = (1.0 * dpi * height / 4) ** 2
    for i in range(n):
        if radius is not None:
            s = (scaler * radius[i]) ** 2
            if s > maxs:
                s = maxs
            plt.scatter(pos[i, 0], pos[i, 1], color='black', s=s, lw=0)
        else:
            plt.scatter(pos[i, 0], pos[i, 1], color='black', s=50, lw=0)

    if label is not None:
        for each in label:
            ax.text(pos[each, 0], pos[each, 1], each, fontsize=5)

    # Plot the edges
    if edges is not None:
        segments = [[pos[i, :], pos[j, :]] for (i, j) in edges]
        lc = LineCollection(segments, colors="lightgray", zorder=0, cmap=plt.cm.Reds)
        lc.set_linewidths(np.full(len(segments), 0.5))
        ax.add_collection(lc)

    if label != None:
        for i in range(n):
            ax.annotate(label[i], (pos[i, 0], pos[i, 1]))

    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    plt.savefig('{}.pdf'.format(name), bbox_inches='tight', pad_inches=0.0)
    plt.clf()


def add_radius(D, nodeweight, scaler=0.03):
    avgsize = np.average(nodeweight)
    avgdist = np.average(D)
    scale = scaler * avgdist / np.sqrt(avgsize)
    n = len(nodeweight)
    r = [scale * np.sqrt(nodeweight[i]) for i in range(n)]  # node radius
    # update Distance by adding radius
    for i in range(n):
        for j in range(n):
            D[i, j] = D[i, j] + r[i] + r[j]

    return D, r

def get_ppr(G, t, alpha=0.5):
    G = G.to_undirected()
    A = nx.adjacency_matrix(G)
    P = preprocessing.normalize(A, norm='l1', axis=1)
    I = identity(P.shape[0])
    ppr = I
    for i in range(t):
        ppr = (1 - alpha) * P.dot(ppr) + I
    ppr = alpha * ppr

    return ppr

def get_pprdegree_dist(G=None, is_vanilla=False, is_vanilla_scale=False, PPR=None, t=15):
    if PPR is None:
        assert (G is not None)
        PPR = get_ppr(G, t=t, alpha=0.2)
        PPR = PPR.todense()
        n = PPR.shape[0]
        if is_vanilla==False:
            for u in range(n):
                PPR[u] = PPR[u] * (G.degree(u) + 1)
        ppr = PPR + PPR.T
    else:
        PPR = PPR.todense()
        n = PPR.shape[0]
        ppr = PPR + PPR.T

    if is_vanilla_scale == False:
        maxval = 2 * np.log2(n)
        minval = 2
        ppr[ppr == 0] = pow(2, 1 - maxval)
        D = 1 - np.log2(ppr)
        cnt = 0
        for u in range(n):
            for v in range(n):
                if (u!=v and (D[u, v]<minval or D[u, v]>maxval)):
                    cnt+=1
                D[u, v] = min(max(minval, D[u, v]), maxval)
        # print ("boundary rate:{:.2f} ".format(1.0*cnt/n/n))
    else:
        maxval = 2 * np.log2(n)
        ppr[ppr == 0] = pow(2, 1 - maxval)
        D = 1 - np.log2(ppr)

    return D


def pprdegree_relativemds(G=None, PPR=None, nodeweight=None, is_vanilla=False, is_vanilla_scale=False, t=15):
    D = get_pprdegree_dist(G, is_vanilla=is_vanilla,is_vanilla_scale=is_vanilla_scale, PPR=PPR, t=t)
    # for nd in leftbottm:
    #     print (nd, D[nd])
    if nodeweight is not None:
        D, r = add_radius(D, nodeweight)
    W = 1 / D / D
    scalar = 2

    # print ("PPRDist mean:{} std:{}".format(np.mean(D),np.std(D)))

    if G is not None:
        n = len(G)
        ratio = 0
        candidate = []
        for c in sorted(nx.connected_components(G), key=len):
            ratio += len(c) / (1.0 * n)
            if ratio <= 0.2:
                candidate.append(c)
            else:
                break
        # print(candidate)
        if len(candidate) > 1:
            for a in candidate:
                for b in candidate:
                    if a == b:
                        continue
                    for (i, j) in [(x, y) for x in a for y in b]:
                        W[i, j] *= scalar

    # mds = relativeMDS(n_components=2,  max_iter=3000, eps=1e-9, random_state=0, dissimilarity='precomputed', n_jobs=1)
    mds = relativeMDS(n_components=2, random_state=0, dissimilarity='precomputed', n_jobs=1)
    pos = mds.fit_transform(D, weights=W)
    # pca = PCA(n_components=2,svd_solver="full")
    # pos = pca.fit_transform(pos)

    if nodeweight is not None:
        return pos, r
    else:
        return pos


def dist_mds(G):
    start = time.time()
    n, m = G.number_of_nodes(), G.number_of_edges()
    len_dict = dict(nx.all_pairs_shortest_path_length(G))

    row_ids = []
    col_ids = []
    dists = []
    for i in range(n):
        for j in range(n):
            row_ids.append(i)
            col_ids.append(j)
            if i == j or j not in len_dict[i]:
                dist = 0
            else:
                dist = len_dict[i][j] + 1
            dists.append(dist)

    D = csr_matrix((dists, (row_ids, col_ids)))
    mds = manifold.MDS(n_components=2, random_state=0, dissimilarity='precomputed', n_jobs=1)
    pos = mds.fit_transform(D.toarray())
    t = time.time() - start

    return t, pos
