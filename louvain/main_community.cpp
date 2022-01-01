#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <ctime>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <numeric>
#include <set>
#include <algorithm>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "graph_binary.h"
#include "community.h"
#define LOUVAIN 0
#define LOUVAINMERGE 1
using namespace std;
namespace bp = boost::property_tree;

char* filename = NULL;
char *filename_w = NULL;
char *filename_part = NULL;
int type = UNWEIGHTED;
int algorithm = -1;
int fid = -1;
int seed = -1;
map<int, string> filelist = {{4,"trust"},{5,"scinet"},{6, "amazon"},{7, "youtube"},
                             {8, "dblp"},{9, "orkut"},{10, "it"},{11, "tw"}};

double precision = 0.000001;
int display_level = -2;
int k = 25;
int thread_num = 1;
int chunk_size = 32;
double total_time = 0;
bool verbose = false;
bool output = false;
string rootpath = "../";

// a series of node partitions; key is community value is node
// structure: component(levels(community(nodes)))
vector<vector<vector<vector<int>>>> partitions;
// small connected component no need louvainl
vector<vector<int>> smallpartition;



double time_by(double &start);
double louvain(Community &c);
vector<vector<vector<int>>>& louvainPlus(Community &c,vector<vector<vector<int>>> &partit);
void write_partition(const string& hiename, const string& mapname, const string &rootname);
vector<vector<int>>& scc(Graph &G,vector<vector<int>> &SCCG);

double time_by(double &start){
    return (omp_get_wtime()-start);
}

int
parseLine(char *line) {
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char *p = line;
    while (*p < '0' || *p > '9') p++;
    line[i - 3] = '\0';
    i = atoi(p);
    return i;
}

void
usage(char *prog_name, const char *more) {
    cerr << more;
    cerr << "usage: " << prog_name
         << "[-f file_no] [-a algorithm] [-k partition_size] [-v] [-o]" << endl << endl;
    cerr << "-f: file containing the graph to decompose in communities." << endl;
    cerr << "-a: 0 is louvain or 1 is louvainPlus." << endl;
    cerr << "-k: threshold of partition size." << endl;
    cerr << "-v: verbose mode." << endl;
    cerr << "-s: random seed." << endl;
    cerr << "-o: output the partition." << endl;
    exit(0);
}

void
parse_args(int argc, char **argv) {
    if (argc < 2)
        usage(argv[0], "Bad arguments number\n");

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'a': // input arg to select different algorithm: louvain or louvainPlus
                    algorithm = atoi(argv[i + 1]);
                    i++;
                    break;
                case 'f': // input arg to select input file id
                    fid = atoi(argv[i + 1]);
                    i++;
                    break;
                case 'k':
                    k = atoi(argv[i + 1]);
                    i++;
                    break;
                case 'o':
                    output = true;
                    i++;
                    break;
                case 'n':
                    thread_num = atoi(argv[i + 1]);
                    i++;
                    break;
                case 'v':
                    verbose = true;
                    i++;
                    break;
                case 's':
                    seed = atoi(argv[i + 1]);
                    i++;
                    break;
                default:
                    usage(argv[0], "Unknown option\n");
            }
        } else {
            usage(argv[0], "Unknown option\n");
        }
    }
}

void
display_time(const char *str) {
    time_t rawtime;
    time(&rawtime);
    cerr << str << ": " << ctime(&rawtime);
}


int
main(int argc, char **argv) {
//    srand(time(NULL) + getpid());
//    srand(1);
    parse_args(argc, argv);
    double starttime;

    string strFname = rootpath+"dataset/"+filelist[fid] +".bin";
    filename = &strFname[0];
    cerr << filename << endl;

    if (algorithm==LOUVAIN){
        Community c(filename, filename_w, type, -1, precision,1);
        starttime = omp_get_wtime();
        louvain(c);
    }
    else if(algorithm==LOUVAINMERGE){
        Graph g(filename, filename_w, type);
        starttime = omp_get_wtime();
        vector<vector<int>> sccs;
        scc(g,sccs);
        cerr << "SCC time: "<< time_by(starttime) << endl;
//        int chunks[] = {1};
//        int threads[] = {64,32,16,8,4,2,1};
        int threads[] = {1};
//        int threads[] = {64};
        for(int trs:threads){
            thread_num = trs<omp_get_max_threads()? trs:omp_get_max_threads();
            cerr<<"chunk: "<<chunk_size<<" threads: "<<thread_num<<endl;
            initmaps(g.nb_nodes,thread_num);
            total_time = 0;
            if(sccs.size()==1){
                double start = omp_get_wtime();
                Community c(g, -1, precision,thread_num,chunk_size);
                total_time+=time_by(start);
                vector<vector<vector<int>>> partit;
                louvainPlus(c,partit);
                partitions.emplace_back(partit);
            } else{
                int maxsize = 0;
                int idx = -1;
                for (int i = 0; i < sccs.size(); ++i) {
                    if(sccs[i].size()>maxsize){
                        maxsize=sccs[i].size();
                        idx = i;
                    }
                }
                Graph subg(g,sccs[idx]);
                double start = omp_get_wtime();
                Community c(subg, -1, precision,thread_num,chunk_size);
                total_time+=time_by(start);
                // the hierarchy structure for the largest component
                vector<vector<vector<int>>> partit;
                louvainPlus(c,partit);
                vector<vector<int>>& coarsest = partit.back();
                if(coarsest.size()>1){
                    vector<int> child(coarsest.size(),0);
                    iota(std::begin(child), std::end(child), 0);
                    vector<vector<int>> root = {child};
                    partit.emplace_back(root);
                }
                partitions.emplace_back(partit);
            }
            cout <<total_time<<endl;
            if (thread_num==1){
                if (output){
                    string hiename = rootpath+"louvain/hierachy-output/"+filelist[fid] +"_"+to_string(k)+".dat";
                    string rootname = rootpath+"louvain/hierachy-output/"+filelist[fid] +"_"+to_string(k)+".root";
                    string mapname = rootpath+"louvain/mapping-output/"+filelist[fid] +"_"+to_string(k)+".dat";
                    write_partition(hiename, mapname, rootname);
                }
            }
            partitions.clear();
            smallpartition.clear();
        }
    }
    else{
        cerr << "error input algorithm: can only select 0(Louvain) or 1(LouvainMerge)." << endl;
        exit(-1);
    }
}

vector<vector<int>>& scc(Graph &G,vector<vector<int>> &SCCG){
    int vnums = G.nb_nodes;
    vector<int> DFN(vnums,0), q(vnums, 0), LOW(vnums,0), st(vnums,0), viter(vnums,0), id(vnums,0);
    vector<bool> visit(vnums, false);
    int top, Tindex, stop, snum=0;
    top=-1; Tindex=0;stop=-1;
    for(int u=0;u<vnums;u++){
        if(!DFN[u]){
            int v;
            DFN[u] = LOW[u] = ++Tindex;
            visit[u] = true;
            q[++top] = u;
            st[++stop] =u;
            while (stop!=-1){
                u = st[stop];
                pair<vector<int>::iterator, vector<float>::iterator > p = G.neighbors(u);
                for(int iter = viter[u];iter < G.nb_neighbors(u);iter++) {
                    v = *(p.first+iter);
                    if (!DFN[v]) {
                        DFN[v] = LOW[v] = ++Tindex;
                        visit[v] = true;
                        q[++top] = v;
                        st[++stop] = v;
                        break;
                    } else if (DFN[v] > DFN[u] && LOW[v] < LOW[u]) {
                        LOW[u] = LOW[v];
                    } else if (visit[v] && DFN[v] < LOW[u])
                        LOW[u] = DFN[v];
                }
                if(u == st[stop]) {
                    if (DFN[u] == LOW[u]) {
                        do {
                            v = q[top--];
                            visit[v] = false;
                            id[v] = snum;
                        } while (v != u);
                        snum++;
                    }
                    stop--;
                }
            }
        }
    }
    SCCG.resize(snum);
    DFN.assign(vnums,0);
    for(int i=0;i<vnums;i++){
        SCCG[id[i]].emplace_back(i);
    }
    return SCCG;
}

vector<vector<vector<int>>>& louvainPlus(Community &c,vector<vector<vector<int>>> &partit){

    Graph g;
    double mod = c.modularity(), new_mod;
    int level = 0;

    double mergetime = 0;
    double graphtime = 0;
    double misctime = 0;
    double start;

    while(c.g.nb_nodes>k) {
        if (verbose) {
            cerr << "level " << level << ":\n";
//            display_time("  start computation");
            cerr << "  network size: "
                 << c.g.nb_nodes << " nodes, "
                 << c.g.nb_links << " links, "
                 << c.g.total_weight << " weight."<< endl;
        }
        start = omp_get_wtime();
        c.one_level_new(k,seed);
        new_mod = c.modularity();
        mergetime+=time_by(start);

        if (verbose)
            cerr << "  modularity increased from " << mod << " to " << new_mod << endl;

        if (++level == display_level)
            g.display();
        if (display_level == -1)
            c.display_partition();

        vector<vector<int>> comm_nodes;
        start = omp_get_wtime();
        g = c.partition2graph_binary(comm_nodes);
        graphtime+=time_by(start);
        int maxsize=0;
        int minsize=100000000;
        int avgsize=0;
        for (auto comm:comm_nodes){
            int s = comm.size();
            if (s<minsize) minsize = s;
            if (s>maxsize) maxsize = s;
            avgsize += s;
        }
        cout<<"community size: min: "<<minsize<<" max: "<<maxsize<<" avg: "<<avgsize/comm_nodes.size()<<endl;
        partit.push_back(comm_nodes);
        start = omp_get_wtime();
        c = Community(g, -1, precision,thread_num,chunk_size);
        mod = new_mod;
        misctime += time_by(start);

    }
    if (verbose) {
        cerr << "level " << level << ":\n";
        cerr << "  network size: "
             << c.g.nb_nodes << " nodes, "
             << c.g.nb_links << " links, "
             << c.g.total_weight << " weight."<< endl;
        cerr <<"merge time: "<<mergetime<<" generate graph time: "<<graphtime<<endl;
    }
    total_time += mergetime+graphtime+misctime;
    return partit;
}
double louvain(Community &c) {
    Graph g;
    bool improvement = true;
    double mod = c.modularity(), new_mod;
    int level = 0;
    do {
        if (verbose) {
            cerr << "level " << level << ":\n";
            display_time("  start computation");
            cerr << "  network size: "
                 << c.g.nb_nodes << " nodes, "
                 << c.g.nb_links << " links, "
                 << c.g.total_weight << " weight." << endl;
        }

        improvement = c.one_level();
        new_mod = c.modularity();
        if (verbose)
            cerr << "  modularity increased from " << mod << " to " << new_mod << endl;
        if (++level == display_level)
            g.display();
        if (display_level == -1)
            c.display_partition();
        g = c.partition2graph_binary();
        c = Community(g, -1, precision,1,1);

        mod = new_mod;
        if (verbose)
            display_time("  end computation");

        if (filename_part != NULL && level == 1) // do at least one more computation if partition is provided
            improvement = true;
    } while (improvement);
    return new_mod;
}


void write_partition(const string& hiename, const string& mapname, const string &rootname){
    std::ofstream mapofs,hierofs,rotofs;
    mapofs.open(mapname);
    hierofs.open(hiename);
    rotofs.open(rootname);
    assert(mapofs.is_open() && hierofs.is_open() && rotofs.is_open());
    int id = 0;
    for(;id<partitions.size();id++){
        vector<vector<vector<int>>>& partit=partitions[id];
        vector<vector<int>>& l1 = partit[0];
        vector<vector<int>> leafs=l1;
        for (int i = 0; i < l1.size(); ++i) {
            string supernode = "c"+to_string(id)+"_l"+to_string(1)+"_"+to_string(i);
            mapofs<<supernode<<endl;
            mapofs<<l1[i].size()<<endl;
            for (int ele:l1[i])
                mapofs<<ele<<endl;
        }
        // for higher level, bottom up store mapping file and hierachy file
        for (int l = 1; l < partit.size(); ++l) {
            vector<vector<int>> &p = partit[l];
            vector<vector<int>> newleaf;
            newleaf.resize(p.size());
            for (int i = 0; i < p.size(); ++i) {
                vector<int> &childsupernode = p[i];
                string supernode = "c"+to_string(id)+"_l" + to_string(l+1) + "_" + to_string(i);
                if(l==partit.size()-1)
                    rotofs<<supernode<<endl;
                hierofs<<supernode<<endl;
                hierofs<<childsupernode.size()<<endl;
                for (int ele:childsupernode)
                    hierofs<<ele<<endl;
                vector<int> curleaf;
                for(int each:childsupernode){
                    curleaf.insert(curleaf.end(),leafs[each].begin(),leafs[each].end());
                }
                mapofs<<supernode<<endl;
                mapofs<<curleaf.size()<<endl;
                for (int ele:curleaf)
                    mapofs<<ele<<endl;
                newleaf[i] = curleaf;
            }
            leafs = newleaf;
        }
        id++;
    }
    for (auto & partit : smallpartition) {
        string supernode = "c"+to_string(id)+"_l1_0";
        rotofs<<supernode<<endl;
        mapofs<<supernode<<endl;
        mapofs<<partit.size()<<endl;
        for (int ele:partit)
            mapofs<<ele<<endl;
        id++;
    }
    mapofs.close();
    hierofs.close();
}

