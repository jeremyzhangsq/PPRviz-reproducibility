#include <fstream>
#include "graph_binary.h"

Graph::Graph() {
    nb_nodes = 0;
    nb_links = 0;
    total_weight = 0;
}


Graph::Graph(Graph& orginG, vector<int> &nodes){
    nb_nodes = nodes.size();
    int oldsize = orginG.nb_nodes;
    // map orginal node id to new node id
    old2new.resize(oldsize, -1);
    // map new node id to orginal node id
    new2old.resize(nb_nodes,0);
    for(int node : nodes)
        old2new[node]=0;
    int final = 0;
    for (int i = 0; i < oldsize; i++){
        if (old2new[i] != -1){
            old2new[i] = final;
            new2old[final] = i;
            final++;
        }
    }
    // shrink the link, degree and weights by nodes
    vector<int> &oldlinks = orginG.links;
    vector<long> &olddegrees = orginG.degrees;
    vector<float> &oldweights = orginG.weights;
    degrees.resize(nb_nodes);
    total_weight=0;
    int curnode = 0;
    int ngbcnt = 0;
    long deg = 0;
    int ngbsize = orginG.nb_neighbors(curnode);
    // if no weights assigned, only shrink link and degrees
    if (oldweights.empty()){
        for (int i = 0; i < orginG.nb_links;) {
            // if current node not belong to subgraph, then skip its neighbors
            int newnode = old2new[curnode];
            if (newnode == -1){
                i += orginG.nb_neighbors(curnode++);
                if(curnode<oldsize)
                    ngbsize = orginG.nb_neighbors(curnode);
                continue;
            }
            int newid = old2new[oldlinks[i]];
            ngbcnt++;
            links.push_back(newid);
            // increase the degree
            deg++;

            if(ngbcnt==ngbsize){
                // update the deg for current node and move to next
                // make the degree accumated
                if (newnode==0)
                    degrees[newnode] = deg;
                else
                    degrees[newnode] = degrees[newnode-1]+deg;
                deg=0,ngbcnt=0;
                if(++curnode<oldsize)
                    ngbsize = orginG.nb_neighbors(curnode);
            }
            i++;
        }
    }
    else{
        for (int i = 0; i < orginG.nb_links;) {
            // if current node not belong to subgraph, then skip its neighbors
            int newnode = old2new[curnode];
            if (newnode == -1){
                i += orginG.nb_neighbors(curnode++);
                if(curnode<oldsize)
                    ngbsize = orginG.nb_neighbors(curnode);
                continue;
            }
            int newid = old2new[oldlinks[i]];
            ngbcnt++;

            links.push_back(newid);
            weights.push_back(oldweights[i]);
            // increase the degree
            deg++;

            if(ngbcnt==ngbsize){
                // update the deg for current node and move to next
                // make the degree accumated
                if (newnode==0)
                    degrees[newnode] = deg;
                else
                    degrees[newnode] = degrees[newnode-1]+deg;
                deg=0,ngbcnt=0;
                if(++curnode<oldsize)
                    ngbsize = orginG.nb_neighbors(curnode);
            }
            i++;
        }
    }

    assert(degrees[nb_nodes-1]==links.size());
    // update nb_links total_weight
    nb_links=links.size();
    for (int i = 0; i < nb_nodes; i++) {
        total_weight += (double) weighted_degree(i);
    }

}

Graph::Graph(char *filename, char *filename_w, int type) {
    ifstream finput;
    finput.open(filename, fstream::in | fstream::binary);

    // Read number of nodes on 4 bytes
    finput.read((char *) &nb_nodes, 4);
    assert(finput.rdstate() == ios::goodbit);

    // Read cumulative degree sequence: 8 bytes for each node
    // cum_degree[0]=degree(0); cum_degree[1]=degree(0)+degree(1), etc.
    degrees.resize(nb_nodes);
    finput.read((char *) &degrees[0], nb_nodes * 8);

    // Read links: 4 bytes for each link (each link is counted twice)
    nb_links = degrees[nb_nodes - 1];
    links.resize(nb_links);
    finput.read((char *) (&links[0]), (long) nb_links * 8);

    // IF WEIGHTED : read weights: 4 bytes for each link (each link is counted twice)
    weights.resize(0);
    total_weight = 0;
    if (type == WEIGHTED) {
        ifstream finput_w;
        finput_w.open(filename_w, fstream::in | fstream::binary);
        weights.resize(nb_links);
        finput_w.read((char *) &weights[0], (long) nb_links * 4);
    }

    // Compute total weight
    for (int i = 0; i < nb_nodes; i++) {
        total_weight += (double) weighted_degree(i);
    }
}

Graph::Graph(int n, int m, double t, int *d, int *l, float *w) {
/*  nb_nodes     = n;
  nb_links     = m;
  total_weight = t;
  degrees      = d;
  links        = l;
  weights      = w;*/
}


void
Graph::display() {
    for (int node = 0; node < nb_nodes; node++) {
        pair<vector<int>::iterator, vector<float>::iterator> p = neighbors(node);
        cout << node << ":";
        for (int i = 0; i < nb_neighbors(node); i++) {
            if (true) {
                if (weights.size() != 0)
                    cout << " (" << *(p.first + i) << " " << *(p.second + i) << ")";
                else
                    cout << " " << *(p.first + i);
            }
        }
        cout << endl;
    }
}

void
Graph::display_reverse() {
    for (int node = 0; node < nb_nodes; node++) {
        pair<vector<int>::iterator, vector<float>::iterator> p = neighbors(node);
        for (int i = 0; i < nb_neighbors(node); i++) {
            if (node > *(p.first + i)) {
                if (weights.size() != 0)
                    cout << *(p.first + i) << " " << node << " " << *(p.second + i) << endl;
                else
                    cout << *(p.first + i) << " " << node << endl;
            }
        }
    }
}


bool
Graph::check_symmetry() {
    int error = 0;
    for (int node = 0; node < nb_nodes; node++) {
        pair<vector<int>::iterator, vector<float>::iterator> p = neighbors(node);
        for (int i = 0; i < nb_neighbors(node); i++) {
            int neigh = *(p.first + i);
            float weight = *(p.second + i);

            pair<vector<int>::iterator, vector<float>::iterator> p_neigh = neighbors(neigh);
            for (int j = 0; j < nb_neighbors(neigh); j++) {
                int neigh_neigh = *(p_neigh.first + j);
                float neigh_weight = *(p_neigh.second + j);

                if (node == neigh_neigh && weight != neigh_weight) {
                    cout << node << " " << neigh << " " << weight << " " << neigh_weight << endl;
                    if (error++ == 10)
                        exit(0);
                }
            }
        }
    }
    return (error == 0);
}


void
Graph::display_binary(char *outfile) {
    ofstream foutput;
    foutput.open(outfile, fstream::out | fstream::binary);

    foutput.write((char *) (&nb_nodes), 4);
    foutput.write((char *) (&degrees[0]), 4 * nb_nodes);
    foutput.write((char *) (&links[0]), 8 * nb_links);
}
