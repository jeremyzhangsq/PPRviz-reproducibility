#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "lib.h"
using namespace std;

bool sortdesc(const pair<double,int> &a,
              const pair<double,int> &b)
{
    return (a.first > b.first);
}

class Graph {
public:

    vector<vector<int>> g;
    vector<int> deg;
    string data_folder;
    vector<int> nodes;
    int n;
    long long m;
    int max_level;
    double dbar;
    double alpha;
    double delta;
    double pfail;
    double epR;
    double tau;
    int k;
    // sum square of degree
    Graph() = default;
    Graph(string &datapath, double a,int cluster_sze) {
        this->data_folder = datapath;
        init_graph();
        dbar = double(2*m) / double(n);
        alpha = a;
        pfail = 1/(double)n;
        double delta_scale =250.0;
        delta = 1.0/delta_scale/log(n);
        k=cluster_sze;
        tau = 1.0/pow(k*n,0.5);
        max_level=-1;
        epR = 0.5;
        if (verbose)
            cout << "init graph n: " << this->n << " m: " << this->m << endl;
    }

    void init_nm() {
        string attribute_file = data_folder + FILECONNECT + "attribute.txt";
        ifstream attr(attribute_file);
        string line1, line2;
        char c;
        while (true) {
            attr >> c;
            if (c == '=') break;
        }
        attr >> n;
        while (true) {
            attr >> c;
            if (c == '=') break;
        }
        attr >> m;
    }

    void init_graph() {
        init_nm();
        nodes.resize(n);
        deg.resize(n);
        iota (begin(nodes), end(nodes), 0);
        g = vector<vector<int>>(n, vector<int>());
        string graph_file = data_folder+".txt";
        FILE *fin = fopen(graph_file.c_str(), "r");
        int t1, t2;
        while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
            assert(t1 < n);
            assert(t2 < n);
            if(t1 == t2) continue;
            g[t1].push_back(t2);
            g[t2].push_back(t1);
        }
        for (int i = 0; i < n; ++i) {
            int d = g[i].size();
            deg[i] = d;
        }
    }


    double get_avg_degree() const {
        return dbar;
    }


};



#endif
