import networkx as nx
import numpy as np



def get_PPR(G, n=10,mode = "pprdeg"):
    PPR = np.zeros((n, n))
    personal = {i: 0 for i in range(n)}
    for i in range(n):
        each = i
        personal[each] = 1
        tmp = nx.pagerank(G, max_iter=100, tol=1.0e-8, alpha=0.8, personalization=personal)
        for j in range(n):
            PPR[i, j] = tmp[j]
        if mode == "pprdeg":
            PPR[i] *= (G.degree(i))
        personal[each] = 0
    return PPR

# e = [(0,1),(2,3),(2,4),(3,4),(4,6),(5,6),(5,7),(5,8),(6,7),(6,8),(7,8)]
e = [(0,1),(0,2),(0,3),(1,2),(1,3),(1,4),(2,3),(2,4),(2,5),(3,5),(4,5),(5,7),(6,9),(7,8)]
G = nx.Graph(e)
n = 10
PPR = get_PPR(G,n,mode="ppr")
PPRDeg = get_PPR(G,n,mode="pprdeg")
ppr = PPR+PPR.T
UPDS = 1-np.log(ppr)
ppr = PPRDeg+PPRDeg.T
PDS = 1-np.log(ppr)
# for u,v in [(6,0),(6,4),(6,8),(5,7),(7,8)]:
for u, v in [(0,8), (2, 0), (6, 9)]:
    print("{:.2f} & {:.2f} & {:.2f}\\\  \hline".format(PPR[u,v],UPDS[u,v],PDS[u,v]))

print("====================")

cluster = [0,1,2]
leafmapping = {0:[0,3],1:[1,2],2:[5]}
G = nx.Graph(e)
# Gone = nx.Graph([(0,1),(1,2),(0,2)])
Gtwo = nx.Graph()
Gtwo.add_weighted_edges_from([(0,1,4),(1,2,1),(0,2,1)])

PPRraw = get_PPR(G)
# PPRone = get_PPR(Gone,3)
PPRtwo = get_PPR(Gtwo,3)
PPRthree = np.zeros((3, 3))
for i in cluster:
    for j in cluster:
        size = len(leafmapping[i])*len(leafmapping[j])
        val = 0
        PDistval = 0
        for u in leafmapping[i]:
            for v in leafmapping[j]:
               val += PPRraw[u,v]
        PPRthree[i,j] = val/size

for u, v in [(1,0), (2, 0), (2, 1)]:
    print("{:.2f} & {:.2f}\\\  \hline".format(PPRtwo[u,v],PPRthree[u,v]))

print("====================")