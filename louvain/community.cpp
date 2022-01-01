#include "community.h"


vector<iMap<float>> maps;
using namespace std;


Community::Community(char *filename, char *filename_w, int type, int nbp, double minm, int nthread) {
    g = Graph(filename, filename_w, type);
    size = g.nb_nodes;
    thread_num = nthread;
    neigh_map.resize(thread_num);
    neigh_last.resize(thread_num);
    for (int i = 0; i < thread_num; ++i) {
        neigh_map[i].initialize(size);
        neigh_last[i] = 0;
    }


    n2c.resize(size);
    c2n.resize(size);
    in.resize(size);
    tot.resize(size);

    int chunk = (size + thread_num - 1) / thread_num;
#pragma omp parallel num_threads(thread_num)
#pragma omp for schedule(static)
        for (int i = 0; i < thread_num; i++) {
            int t = omp_get_thread_num();
            int start = t * chunk;
            int end = (start + chunk <= size) ? (start + chunk) : size;
            // each thread processes chunk number of supernode
            iMap<float> &map = maps[t];
            for (int node = start; node < end; node++) {
                n2c[node] = node;
                c2n[node] = {node};
                tot[node] = g.weighted_degree(node);
                in[node] = g.nb_selfloops(node);
            }
        }

    nb_pass = nbp;
    min_modularity = minm;
}


Community::Community(Graph gc, int nbp, double minm, int nthread, int chunk) {
    g = gc;
    size = g.nb_nodes;
    thread_num = nthread;
    chunksize = chunk;
    neigh_map.resize(thread_num);
    neigh_last.resize(thread_num);
    for (int i = 0; i < thread_num; ++i) {
        neigh_map[i].initialize(size);
        neigh_last[i] = 0;
    }

    n2c.resize(size);
    c2n.resize(size);
    in.resize(size);
    tot.resize(size);

    int chunks = (size + thread_num - 1) / thread_num;
#pragma omp parallel num_threads(thread_num)
#pragma omp for schedule(static)
    for (int i = 0; i < thread_num; i++) {
        int t = omp_get_thread_num();
        int start = t * chunks;
        int end = (start + chunks <= size) ? (start + chunks) : size;
        // each thread processes chunk number of supernode
        iMap<float> &map = maps[t];
        for (int node = start; node < end; node++) {
            n2c[node] = node;
            c2n[node] = {node};
            tot[node] = g.weighted_degree(node);
            in[node] = g.nb_selfloops(node);
        }
    }

    nb_pass = nbp;
    min_modularity = minm;
}
void initmaps(int size,int thread_num){
    maps.resize(thread_num);
    for (int i = 0; i < thread_num; ++i) {
        maps[i].initialize(size);
    }
}



double
Community::modularity() {
    double q = 0.;
    double m2 = (double) g.total_weight;

    for (int i = 0; i < size; i++) {
        if (tot[i] > 0)
            q += (double) in[i] / m2 - ((double) tot[i] / m2) * ((double) tot[i] / m2);
    }

    return q;
}

void
Community::neigh_comm(int node) {
    neigh_last[0] = 0;
    neigh_map[0].clean();
    pair<vector<int>::iterator, vector<float>::iterator> p = g.neighbors(node);
    int deg = g.nb_neighbors(node);
    iMap<double> &nmap = neigh_map[0];
    nmap.insert(n2c[node],0);
    neigh_last[0] = 1;

    for (int i = 0; i < deg; i++) {
        int neigh = *(p.first + i);
        int neigh_comm = n2c[neigh];
        double neigh_w = (g.weights.size() == 0) ? 1. : *(p.second + i);

        if (neigh != node) {
            if (nmap.notexist(neigh_comm)) {
                nmap.insert(neigh_comm,0);
                neigh_last[0]++;
            }
            nmap[neigh_comm] += neigh_w;
        }
    }
}
void
Community::neigh_comm(int comm_id, vector<int> &community, int tid) {
    int &last = neigh_last[tid];
    iMap<double> &nmap = neigh_map[tid];
    nmap.clean();
    last = 0;
    nmap.insert(comm_id,0);
    last = 1;
    for(int node:community){
        pair<vector<int>::iterator, vector<float>::iterator> p = g.neighbors(node);
        int deg = g.nb_neighbors(node);

        for (int i = 0; i < deg; i++) {
            int neigh = *(p.first + i);
            int neigh_comm = n2c[neigh];
            assert(neigh_comm!=-1);
            double neigh_w = (g.weights.size() == 0) ? 1. : *(p.second + i);

            if (nmap.notexist(neigh_comm)) {
                nmap.insert(neigh_comm,0);
                last++;
            }
            nmap[neigh_comm] += neigh_w;
        }
    }
}


void
Community::display_partition() {
    vector<int> renumber(size, -1);
    for (int node = 0; node < size; node++) {
        renumber[n2c[node]]++;
    }

    int final = 0;
    for (int i = 0; i < size; i++)
        if (renumber[i] != -1)
            renumber[i] = final++;

    for (int i = 0; i < size; i++)
        cout << i << " " << renumber[n2c[i]] << endl;
}

Graph
Community::partition2graph_binary() {
    // Renumber communities
    vector<int> renumber(size, -1);
    for (int node=0 ; node<size ; node++) {
        renumber[n2c[node]]++;
    }

    int final=0;
    for (int i=0 ; i<size ; i++)
        if (renumber[i]!=-1)
            renumber[i]=final++;

    // Compute communities
    vector<vector<int> > comm_nodes(final);
    for (int node=0 ; node<size ; node++) {
        comm_nodes[renumber[n2c[node]]].push_back(node);
    }

    // Compute weighted graph
    Graph g2;
    g2.nb_nodes = comm_nodes.size();
    g2.degrees.resize(comm_nodes.size());

    int comm_deg = comm_nodes.size();
    for (int comm=0 ; comm<comm_deg ; comm++) {
        map<int,float> m;
        map<int,float>::iterator it;

        int comm_size = comm_nodes[comm].size();
        for (int node=0 ; node<comm_size ; node++) {
            pair<vector<int>::iterator, vector<float>::iterator> p = g.neighbors(comm_nodes[comm][node]);
            int deg = g.nb_neighbors(comm_nodes[comm][node]);
            for (int i=0 ; i<deg ; i++) {
                int neigh        = *(p.first+i);
                int neigh_comm   = renumber[n2c[neigh]];
                double neigh_weight = (g.weights.size()==0)?1.:*(p.second+i);

                it = m.find(neigh_comm);
                if (it==m.end())
                    m.insert(make_pair(neigh_comm, neigh_weight));
                else
                    it->second+=neigh_weight;
            }
        }
        g2.degrees[comm]=(comm==0)?m.size():g2.degrees[comm-1]+m.size();
        g2.nb_links+=m.size();


        for (it = m.begin() ; it!=m.end() ; it++) {
            g2.total_weight  += it->second;
            g2.links.push_back(it->first);
            g2.weights.push_back(it->second);
        }
    }

    return g2;
}


Graph
Community::partition2graph_binary(vector<vector<int>> &comm_nodes) {
    // Renumber communities
    vector<int> renumber(size, -1);
    for (int node = 0; node < size; node++) {
        renumber[n2c[node]]++;
    }

    int final = 0;
    for (int i = 0; i < size; i++)
        if (renumber[i] != -1)
            renumber[i] = final++;

    // Compute communities
    comm_nodes.resize(final);
    for (int node = 0; node < size; node++) {
        comm_nodes[renumber[n2c[node]]].push_back(node);
    }

    // Compute weighted graph
    Graph g2;
    g2.nb_nodes = comm_nodes.size();
    g2.degrees.resize(comm_nodes.size());
    int comm_deg = comm_nodes.size();
    vector<int> degrees(comm_deg,0);
    vector<float> weights(comm_deg,0);
    int chunk = (comm_deg + thread_num - 1) / thread_num;
    //printf("%llu %llu %d %d\n ", n, m, p, num_chunk);
#pragma omp parallel num_threads(thread_num)
    {
        vector<float> weights_private;
        vector<int> links_private;
#pragma omp for schedule(static)
        for (int i = 0; i < thread_num; i++) {
            int t = omp_get_thread_num();
            int start = t * chunk;
            int end = (start + chunk <= comm_deg) ? (start + chunk) : comm_deg;
            // each thread processes chunk number of supernode
            iMap<float> &map = maps[t];
            for (int comm = start; comm < end; comm++) {
                int comm_size = comm_nodes[comm].size();
                for (int node = 0; node < comm_size; node++) {
                    pair<vector<int>::iterator, vector<float>::iterator> p = g.neighbors(comm_nodes[comm][node]);
                    int deg = g.nb_neighbors(comm_nodes[comm][node]);
                    for (int i = 0; i < deg; i++) {
                        int neigh = *(p.first + i);
                        int neigh_comm = renumber[n2c[neigh]];
                        double neigh_weight = (g.weights.size() == 0) ? 1. : *(p.second + i);

                        if (!map.exist(neigh_comm))
                            map.insert(neigh_comm, neigh_weight);
                        else
                            map[neigh_comm] += neigh_weight;
                    }
                }
                long msize = map.occur.m_num;
                degrees[comm] = msize;
                float com_weight = 0;
                for (int i = 0; i < msize; i++) {
                    int node = map.occur[i];
                    float weight = map[node];
                    links_private.push_back(node);
                    weights_private.push_back(weight);
                    com_weight+= weight;
                }
                weights[comm] = com_weight;
                map.clean();
            }
        }
#pragma omp for schedule(static) ordered
        for (int i = 0; i < thread_num; i++) {
#pragma omp ordered
            {
                g2.weights.insert(g2.weights.end(), weights_private.begin(), weights_private.end());
                g2.links.insert(g2.links.end(), links_private.begin(), links_private.end());
            }
        }
    }
    for (int comm = 0; comm < comm_deg; ++comm) {
        g2.degrees[comm] = (comm == 0) ? degrees[comm] : g2.degrees[comm - 1] + degrees[comm];
        g2.nb_links += degrees[comm];
        g2.total_weight += weights[comm];

    }
    // map back to old id for level 1
    if(!g.new2old.empty()){
        for(vector<int> &list:comm_nodes){
            for (int i = 0; i < list.size(); ++i) {
                list[i] = g.new2old[list[i]];
            }
        }
    }

    return g2;
}


void
Community::one_level_new(const int &k, int seed) {
    vector<int> Si(g.nb_nodes);
    std::iota(std::begin(Si), std::end(Si), 0);
    int left = 0;
//    random_device rd;
//    mt19937 rng(rd());
    mt19937 rng(seed);
    shuffle ( Si.begin(), Si.end(),rng );
    vector<int> inqueue(g.nb_nodes,true);
    // initial variable for each thread
    vector<vector<int>> candidates(thread_num);
    vector<vector<double>> w_degrees(thread_num);
    vector<vector<int>> best_neighbors(thread_num);
    vector<vector<int>> best_links(thread_num);

    for (int i = 0; i < thread_num; ++i) {
        candidates[i].reserve(chunksize);
        w_degrees[i].reserve(chunksize);
        best_neighbors[i].reserve(chunksize);
        best_links[i].reserve(chunksize);
    }
    vector<int> exist_neighbor(g.nb_nodes,0);
    vector<int> invalid(g.nb_nodes,false);
    vector<int> invalid_name; invalid_name.reserve(chunksize*thread_num*2);
    int cnt = 0;
    while (left < (int) Si.size()) {
        for (int j = 0; j < chunksize; ++j) {
            int i = 0;
            while ( i < thread_num) {
                if (left == (int) Si.size()) break;
                if (!inqueue[Si[left]]){
                    left++;
                    i++;
                    continue;
                }
                candidates[i].push_back(Si[left]);
                left++;
                i++;
            }
            if (left == (int) Si.size()) break;
        }

        cnt ++;
//        double begin = omp_get_wtime();
#pragma omp parallel num_threads(thread_num)
#pragma omp for schedule(static)
        for (int i = 0; i < thread_num; ++i) {
            int id = omp_get_thread_num();
            for(int v_comm:candidates[id]){
//                assert(v_comm!=-1);
                vector<int> &comm = c2n[v_comm];
                int v_comm_size = comm.size();
                // the sum of the weightes of the links incident to current community
                double w_degree = g.weighted_degree(comm);
                // computation of all neighboring communities of current community
                neigh_comm(v_comm, comm, id);
                iMap<double> &nmap = neigh_map[id];
                if (nmap.occur.m_num == 1) {
                    continue;
                }

                int &last = neigh_last[id];

                double node_nblink = nmap[v_comm];
//                assert(tot[v_comm] == w_degree);
//                assert(in[v_comm] == node_nblink);
                int u_comm = nmap.occur[1];
                int u_comm_size = c2n[u_comm].size();
                double best_nblinks = nmap[u_comm];
                double best_increase = modularity_gain(u_comm, best_nblinks, w_degree);;

                for (int i = 0; i < last; i++) {
                    int w_comm = nmap.occur[i];
                    if (v_comm == w_comm)
                        continue;
                    int w_comm_size = c2n[w_comm].size();
                    double increase = modularity_gain(w_comm, nmap[w_comm], w_degree);
                    if (v_comm_size + w_comm_size <= k) {
                        if (increase > best_increase || v_comm_size + u_comm_size > k) {
                            u_comm = w_comm;
                            u_comm_size = w_comm_size;
                            best_nblinks = nmap[w_comm];
                            best_increase = increase;
                        }
                    }
                    else if (w_comm_size < u_comm_size) {
                        u_comm = w_comm;
                        u_comm_size = w_comm_size;
                        best_nblinks = nmap[w_comm];
                        best_increase = increase;
                    }
                }
                best_neighbors[id].push_back(u_comm);
                best_links[id].push_back(best_nblinks);
                w_degrees[id].push_back(w_degree);
            }
        }

        // update tot and in for ngbr_comm and remove/insert nodes from/to S
        for(int tid=0;tid<thread_num;tid++){
            vector<int> &bestngbr = best_neighbors[tid];
            if (bestngbr.empty())
                continue;
            vector<int> &candi = candidates[tid];
            vector<double> &wdegree = w_degrees[tid];
            vector<int> &bestnlink = best_links[tid];
//            assert(candi.size()==wdegree.size());
            for(int i=0;i<candi.size();++i){
                int u_comm = bestngbr[i]; //neighbpr
                int v_comm = candi[i];//node itself
                inqueue[v_comm] = false;
                if (invalid[u_comm] || invalid[v_comm]){
                    Si.push_back(v_comm);
                    inqueue[v_comm] = true;
                    continue;
                }
                invalid[u_comm]=true; invalid[v_comm]=true;
                invalid_name.push_back(u_comm);
                invalid_name.push_back(v_comm);
                double real_degree = g.weighted_degree(c2n[v_comm]);
                double w_degree = wdegree[i];
//                assert(real_degree==w_degree);
                int best_nblinks = bestnlink[i];
//                assert(v_comm!=u_comm);
                tot[u_comm] += w_degree;
                in[u_comm] += in[v_comm] + 2 * best_nblinks;
                for(int each:c2n[v_comm]){
                    n2c[each] = u_comm;
                    c2n[u_comm].push_back(each);
                }

                c2n[v_comm].clear();
                tot[v_comm] = 0;
                in[v_comm] = 0;
                if(c2n[u_comm].size()>=k)
                    inqueue[u_comm] = false;
                else{
                    if (inqueue[u_comm] ==false){
                        Si.push_back(u_comm);
                        inqueue[u_comm] = true;
                    }
                }
            }
        }
        for (int each:invalid_name)
            invalid[each]=false;
        invalid_name.clear();
        // empty the candidate list
        for (int i = 0; i < thread_num; ++i) {
            candidates[i].clear();
            candidates[i].reserve(chunksize);
            // clean the data structure
            w_degrees[i].clear();
            candidates[i].reserve(chunksize);
            best_neighbors[i].clear();
            best_links[i].clear();
            candidates[i].reserve(chunksize);
        }
    }

//   cout<<"process iteration: "<<cnt<<endl;
}


bool
Community::one_level() {
    bool improvement=false ;
    int nb_moves;
    int nb_pass_done = 0;
    double new_mod   = modularity();
    double cur_mod   = new_mod;

    vector<int> random_order(size);
    for (int i=0 ; i<size ; i++)
        random_order[i]=i;
    for (int i=0 ; i<size-1 ; i++) {
        int rand_pos = rand()%(size-i)+i;
        int tmp      = random_order[i];
        random_order[i] = random_order[rand_pos];
        random_order[rand_pos] = tmp;
    }

    // repeat while
    //   there is an improvement of modularity
    //   or there is an improvement of modularity greater than a given epsilon
    //   or a predefined number of pass have been done
    do {
        cur_mod = new_mod;
        nb_moves = 0;
        nb_pass_done++;

        // for each node: remove the node from its community and insert it in the best community
        for (int node_tmp=0 ; node_tmp<size ; node_tmp++) {
//      int node = node_tmp;
            int node = random_order[node_tmp];
            int node_comm     = n2c[node];
            double w_degree = g.weighted_degree(node);

            // computation of all neighboring communities of current node
            neigh_comm(node);
            iMap<double> &nmap = neigh_map[0];
            // remove node from its current community
            remove(node, node_comm, nmap[node_comm]);

            // compute the nearest community for node
            // default choice for future insertion is the former community
            int best_comm        = node_comm;
            double best_nblinks  = 0.;
            double best_increase = 0.;
            for (unsigned int i=0 ; i<neigh_last[0] ; i++) {
                double cur_comm = nmap.occur[i];
                double increase = modularity_gain(node,cur_comm, nmap[cur_comm], w_degree);
                if (increase>best_increase) {
                    best_comm     = cur_comm;
                    best_nblinks  = nmap[cur_comm];
                    best_increase = increase;
                }
            }

            // insert node in the nearest community
            insert(node, best_comm, best_nblinks);

            if (best_comm!=node_comm)
                nb_moves++;
        }

        double total_tot=0;
        double total_in=0;
        for (unsigned int i=0 ; i<tot.size() ;i++) {
            total_tot+=tot[i];
            total_in+=in[i];
        }

        new_mod = modularity();
        if (nb_moves>0)
            improvement=true;

    } while (nb_moves>0 && new_mod-cur_mod>min_modularity);

    return improvement;
}

