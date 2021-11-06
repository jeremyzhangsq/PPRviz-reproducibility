import numpy as np
import math
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

# Given three colinear points p, q, r, the function checks if
# point q lies on line segment 'pr'
def onSegment(p, q, r):
    if ((q[0] <= max(p[0], r[0])) and (q[0] >= min(p[0], r[0])) and
            (q[1] <= max(p[1], r[1])) and (q[1] >= min(p[1], r[1]))):
        return True
    return False


def orientation(p, q, r):
    # to find the orientation of an ordered triplet (p,q,r)
    # function returns the following values:
    # 0 : Colinear points
    # 1 : Clockwise points
    # 2 : Counterclockwise

    # See https://www.geeksforgeeks.org/orientation-3-ordered-points/amp/
    # for details of below formula.

    val = (float(q[1] - p[1]) * (r[0] - q[0])) - (float(q[0] - p[0]) * (r[1] - q[1]))
    if (val > 0):

        # Clockwise orientation
        return 1
    elif (val < 0):

        # Counterclockwise orientation
        return 2
    else:

        # Colinear orientation
        return 0


# The main function that returns true if
# the line segment 'p1q1' and 'p2q2' intersect.
def isIntersect(p1, q1, p2, q2):
    # Find the 4 orientations required for
    # the general and special cases
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # General case
    if ((o1 != o2) and (o3 != o4)):
        return True

    # Special Cases

    # p1 , q1 and p2 are colinear and p2 lies on segment p1q1
    if ((o1 == 0) and onSegment(p1, p2, q1)):
        return True

    # p1 , q1 and q2 are colinear and q2 lies on segment p1q1
    if ((o2 == 0) and onSegment(p1, q2, q1)):
        return True

    # p2 , q2 and p1 are colinear and p1 lies on segment p2q2
    if ((o3 == 0) and onSegment(p2, p1, q2)):
        return True

    # p2 , q2 and q1 are colinear and q1 lies on segment p2q2
    if ((o4 == 0) and onSegment(p2, q1, q2)):
        return True

    # If none of the cases
    return False

def get_crossing(G, X):
    res = 0
    edges = G.edges()
    m = len(edges)
    for i in range(m):
        e1 = edges[i]
        for j in range(m):
            if j <= i:
                continue
            e2 = edges[j]
            if e1[0] == e2[0] or e1[1] == e2[1] or e1[0] == e2[1] or e1[1] == e2[0]:
                continue
            if isIntersect(X[e1[0]], X[e1[1]], X[e2[0]], X[e2[1]]):
                res += 1
    return res

# See Section 3.1 in http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.9.3879&rep=rep1&type=pdf
# See Section 3.2 in https://kar.kent.ac.uk/14297/1/graphicalDesignTechniques.pdf
# See https://en.wikipedia.org/wiki/Graph_drawing#cite_note-11
def get_aspect_ratio(X):
    width = X[:, 0].max() - X[:, 0].min()
    height = X[:, 1].max() - X[:, 1].min()
    aspect_graph = width / height
    aspect_view = 4.0 / 3.0
    metric_aspec = max(aspect_graph, aspect_view) / min(aspect_graph, aspect_view) - 1
    return metric_aspec

def cal_angle(vector1, vector2):
    n1 = np.linalg.norm(vector1)
    n2 = np.linalg.norm(vector2)
    if n1 == 0 or n2 == 0:
        return 0
    unit_vector_1 = vector1 / n1
    unit_vector_2 = vector2 / n2
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    if dot_product>=1.0:
        return 0
    angle = np.arccos(dot_product)
    angle = math.fabs(math.degrees(angle))
    return angle

# see Angular Resolution in Section 3.2 in https://kar.kent.ac.uk/14297/1/graphicalDesignTechniques.pdf
def get_angular_resolution(G, X):
    n = X.shape[0]
    metric_angular = 0
    for u in range(n):
        nbrs = G.neighbors(u)
        for vi in nbrs:
            for vj in nbrs:
                if vi >= vj:
                    continue
                pi = X[vi] - X[u]
                pj = X[vj] - X[u]
                angle = cal_angle(pi, pj)
                metric_ij = 1.0 - min(15.0, angle) ** 2 / 15.0 / 15.0
                metric_angular = metric_angular + metric_ij
    return metric_angular

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

def rescale_layout(pos, scale=1.):
    # rescale to [0, scale) in each axis

    # Find max length over all dimensions
    maxlim=0
    for i in range(pos.shape[1]):
        pos[:,i] -= pos[:,i].min() # shift min to zero
        maxlim = max(maxlim, pos[:,i].max())
    if maxlim > 0:
        for i in range(pos.shape[1]):
            pos[:,i] *= scale / maxlim
    return pos


def eva(G, X):
    edges = G.edges()
    try:
        X = rescale_layout(X)
        normX = best_fit_layout(X)
    except: normX = X
    #
    try: a = get_node_distribution(normX)
    except: a = np.infty
    b = uniform_edge_length_coefficient_variance(edges, normX)
    # return a, b
    # add more metrics more revision
    cr = get_crossing(G, normX)
    ar = get_aspect_ratio(normX)
    angr = get_angular_resolution(G, normX)
    return a,b, cr, ar, angr


