

#ifndef INTERACT_FORA_ALGO_H
#define INTERACT_FORA_ALGO_H

#include "lib.h"
#include <boost/serialization/serialization.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
using namespace boost;
using namespace std;


string hiename,mapname,rootname,storepath,rwpath,prpath;
//vector<Bwdidx> offline_idx;
//vector<vector<pair<int,double>>> ppr_idx;
vector<unsigned long long> rw_idx_info_offset;
vector<unsigned long> rw_idx_info_size;
vector<int> bwd_idx_info_offset;
unordered_map<int,pair<int,int>> bwd_idx_info;
vector<float> bwd_idx_in_cluster;
vector<string> root;
unordered_map<string,vector<int>> super2leaf;
vector<vector<int>> level1cluster;
vector<string> hubcluster;
unordered_map<string,vector<int>> super2super;
vector<int> rw_idx;
vector<Fwdidx> fwd_idx;
Bwdidx bwd_idx;
vector<Fwdidx> virtual_fwd_idx;
//vector<Fwdcache> fwd_stack;
vector<vector<int>> leaf2id_stack;
int fwd_stack_top;
int embed_on;
//Fwdcache* fwdcache;
vector<int>* leaf2id;
// current rsum for each leaf source node
//vector<double> rsums;
// store superPPRDeg between supernodes in current level
vector<double> PPRDeg;
Graph graph;
vector<vector<int>> queues;
//string visual_mode; // decide which mode to use: interactive or full?
int thread_num;
double timeElasped;
//vector<double> levelTime;
//vector<double> fwdTime;
//vector<double> rwTime;
vector<float> pr;
string alg;
int isRWIdx = 0;
int isPowerIter = 0;
int isFORASN = 0;
int isFPSN = 0;
int isBPSN = 0;
long long used_rw_size=0;
long long idx_rw_size=0;


inline static double drand(unsigned int &seed){
    return rand_r(&seed)*1.0f/RAND_MAX;
    // return sfmt_genrand_real1(&sfmtSeed);
}

inline static unsigned long lrand(unsigned int &seed) {
    return rand_r(&seed);
    // return sfmt_genrand_uint32(&sfmtSeed);
}

inline int random_walk(int start, unsigned int &seed){
    int cur = start;
    unsigned long k;
    if(graph.g[start].empty()){
        return start;
    }
    while (true) {
        if (drand(seed) < graph.alpha) {
            return cur;
        }
        if (!graph.g[cur].empty()){
            k = lrand(seed)%graph.g[cur].size();
            cur = graph.g[cur][k];
        }
        else{
            cur = start;
        }
    }
}


int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int extract_supernode(string supernode, int &component, int &level){
    size_t pos = 0;
    string token;
    string delimiter = "_";
    int i=1;
    int sid;
    while ((pos = supernode.find(delimiter)) != std::string::npos) {
        token = supernode.substr(1, pos-1);
//        std::cout << token << std::endl;
        if(i==1)
            component = stoi(token);
        else if(i==2)
            level = stoi(token);
        supernode.erase(0, pos + delimiter.length());
        i++;
    }
    assert(i==3);
//    std::cout << supernode << std::endl;
    sid = stoi(supernode);
    return sid;
}

void parse_children(const string &supernode, int &level,vector<string>&partition){
    int component;
    extract_supernode(supernode,component,level);
    if (level<=1)
        return;
    vector<int> &child = super2super[supernode];
    partition.resize(child.size());
    for (int i = 0; i < child.size(); ++i) {
        partition[i] = "c"+to_string(component)+"_l"+to_string(level-1)+"_"+to_string(child[i]);
    }
}

double getMemory(){ //Note: this value is in MB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmRSS:", 6) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result/1024.0;
}
void generate_random_path(vector<string> &path, int maxlevel){

    assert(maxlevel>1);
    int level = maxlevel;
    string item = "c0_l"+to_string(maxlevel)+"_0";
    cerr<<item<<" ";
    path.push_back(item);
    while (level>1){
        int size = super2super[item].size();
        level--;
        item = "c0_l"+to_string(level)+"_"+to_string(super2super[item][rand()%size]);
        cerr<<item<<" ";
        path.push_back(item);
    }
    cerr<<endl;
}
void top_k_hub_cluster(int k){
    // leaf->cluster mapping
    vector<int> leaf2cluster(graph.n,0);
    for (int i = 0; i < level1cluster.size(); ++i) {
        for(auto leaf:level1cluster[i])
            leaf2cluster[leaf]=i;
    }
    // get top-k dnpr node
    priority_queue<pair<float, int>> q;
    hubcluster.reserve(k);
    for (int i = 0; i < pr.size(); ++i) {
        q.push(std::pair<double, int>(pr[i], i));
    }
    for (int i = 0; i < k; ++i) {
        int leaf = q.top().second;
        hubcluster.push_back("c0_l1_"+to_string(leaf2cluster[leaf]));
        q.pop();
    }
//    exit(-1);
}
int load_multilevel(){
    std::ifstream mapifs,hierifs,rotifs;
    mapifs.open(mapname);
    hierifs.open(hiename);
    rotifs.open(rootname);
    string name;
    int max_level=0;
    int size;
    int node;
    int _;
    level1cluster.reserve(graph.n/graph.k);
    assert(mapifs.is_open() && hierifs.is_open() && rotifs.is_open());
    while(!mapifs.eof())
    {
        mapifs >> name >> size;
        vector<int> leafs(size);
        for (int i = 0; i < size; ++i) {
            mapifs >> node;
            leafs[i]=node;
        }
        super2leaf.insert({name,leafs});
        int c,l;
        int id = extract_supernode(name,c,l);
        // the data is stored in order
        if (l==1){
            assert(id==level1cluster.size());
            level1cluster.push_back(leafs);
        }
    }
    while(!rotifs.eof())
    {
        rotifs >> name;
        root.push_back(name);
    }
    while(!hierifs.eof())
    {
        hierifs >> name >> size;
        vector<int> leafs(size);
        for (int i = 0; i < size; ++i) {
            hierifs >> node;
            leafs[i]=node;
        }
        int level;
        extract_supernode(name,_,level);
        if (level>max_level)
            max_level=level;
        super2super.insert({name,leafs});
    }
    mapifs.close();
    hierifs.close();
    if (verbose)
        cout<<"load multilevel success: super2leaf "<<super2leaf.size()<<" super2super "<<super2super.size()<<endl;
    return max_level;
}

inline void serialize_idx(string &path){
    string file_name = path+"_"+alg+".idx";
    std::ofstream ofs(file_name);
    boost::archive::binary_oarchive oa(ofs);
    oa << rw_idx;

    string info_name = path+"_"+alg+".info.offset";
    std::ofstream info_ofs(info_name);
    boost::archive::binary_oarchive info_oa(info_ofs);
    info_oa << rw_idx_info_offset;

    string size_name = path+"_"+alg+".info.size";
    std::ofstream size_ofs(size_name);
    boost::archive::binary_oarchive size_oa(size_ofs);
    size_oa << rw_idx_info_size;
}

void build_rwidx(){
    double start = omp_get_wtime();
    double epsilon =1-exp(-2*graph.epR);
    double Delta;
    if (isFORASN) Delta = graph.k*graph.dbar; // Delta = k*dbar k=25
    else Delta = graph.dbar;
    double deltap = graph.dbar*graph.delta;
//    double deltap = graph.dbar/graph.n;
    double rmax = epsilon*sqrt(deltap*Delta/(2+2*epsilon/3)/2/(double)graph.m/log(1/graph.pfail));
//    rmax *= graph.rmax_scale;
    double omega = (2+2*epsilon/3)*log(1/graph.pfail)/deltap/epsilon/epsilon;
    unsigned int seed = rand();

    // rw_idx = RwIdx( graph.n, vector<int>() );
    rw_idx_info_offset.resize(graph.n);
    rw_idx_info_size.resize(graph.n);
    unsigned long long rw_max_size = 2*graph.m*rmax*omega;
    rw_idx.reserve(rw_max_size);

    {
        unsigned long num_rw;
        for(int source=0; source<graph.n; source++){ //from each node, do rand-walks
            num_rw = ceil(graph.deg[source]*rmax*omega);
            rw_idx_info_offset[source] = rw_idx.size();
            rw_idx_info_size[source] = num_rw;
            for(unsigned long i=0; i<num_rw; i++){ //for each node, do some rand-walks
                int destination = random_walk(source,seed);
                // rw_idx[source].push_back(destination);
                rw_idx.push_back(destination);
            }
        }
    }

    cout <<omp_get_wtime()-start<<endl;
    {
        serialize_idx(rwpath);
    }
#ifdef linux
    cerr << "Memory usage (MB):" << getMemory()<< endl << endl;
#endif
}
void build_rwidx_parallel(){
    double stime = omp_get_wtime();
    double epsilon =1-exp(-2*graph.epR);
    double Delta;
    if (isFORASN) Delta = graph.k*graph.dbar; // Delta = k*dbar k=25
    else Delta = graph.dbar;
    // todo: currently, use dbar*delta as deltap and tau = 1 for randwalk index construction
    double deltap = graph.dbar*graph.delta;
//    double deltap = graph.dbar/graph.n;
    double rmax = epsilon*sqrt(deltap*Delta/(2+2*epsilon/3)/2/(double)graph.m/log(1/graph.pfail));
//    rmax *= graph.rmax_scale;
    double omega = (2+2*epsilon/3)*log(1/graph.pfail)/deltap/epsilon/epsilon;
    unsigned int seed = rand();

    // rw_idx = RwIdx( graph.n, vector<int>() );
    rw_idx_info_offset.resize(graph.n);
    rw_idx_info_size.resize(graph.n);
    unsigned long long rw_max_size = 2*graph.m*rmax*omega;
    rw_idx.reserve(rw_max_size);
//    cout<<rw_max_size<<endl;
    int chunk = (graph.n + thread_num - 1) / thread_num;
//    vector<double> thread_time(thread_num);
    // rw_idx/size/offset for each thread
    vector<vector<int>> rw_idx_private(thread_num);
    vector<vector<int>> rw_idx_size_private(thread_num);
    vector<vector<int>> rw_idx_offset_private(thread_num);
    for (int i = 0; i < thread_num; ++i) {
        rw_idx_private[i].reserve((rw_max_size + thread_num - 1) / thread_num);
        rw_idx_size_private[i].reserve(chunk);
        rw_idx_offset_private[i].reserve(chunk);
    }
    // process list
    vector<int> nodes(graph.n);
    std::iota(std::begin(nodes), std::end(nodes), 0);
    mt19937 rng(1);
    shuffle ( nodes.begin(), nodes.end(),rng );


#pragma omp parallel num_threads(thread_num)
    {
#pragma omp for schedule(static)
        for (int i = 0; i < thread_num; i++) {
            double startt = omp_get_wtime();
            int t = omp_get_thread_num();
            int start = t * chunk;
            int end = (start + chunk <= graph.n) ? (start + chunk) : graph.n;
            for (int i = start; i < end; i++) {//from each node, do rand-walks
                int source = nodes[i];
                unsigned long num_rw;
                num_rw = ceil(graph.deg[source]*rmax*omega);
                rw_idx_size_private[t].push_back(num_rw);
                rw_idx_offset_private[t].push_back(rw_idx_private[t].size());
                for(unsigned long i=0; i<num_rw; i++){ //for each node, do some rand-walks
                    int destination = random_walk(source,seed);
                    // rw_idx[source].push_back(destination);
                    rw_idx_private[t].push_back(destination);
                }
            }
//            thread_time[omp_get_thread_num()] = omp_get_wtime()-startt;
        }
    }

    cout<<omp_get_wtime()-stime<<endl;

//    for (int i = 0; i < thread_num; ++i) {
//        cout <<"thread "<<i<<" rwindextime:"<<thread_time[i]<<endl;
//    }

    // merging
    // generate maps for source and its size and offset
    vector<vector<int>> locations(graph.n); // int[3]: thread-id offset size
    for (int t = 0; t < thread_num; t++) {

        int start = t * chunk;
        int end = (start + chunk <= graph.n) ? (start + chunk) : graph.n;
        for (int i = start; i < end; i++) {//from each node, do rand-walks
            int source = nodes[i];
            locations[source].resize(3);
            locations[source][0] = t;
            locations[source][1] = rw_idx_offset_private[t][i-start];
            locations[source][2] = rw_idx_size_private[t][i-start];
        }
    }
    for (int i = 0; i < graph.n; ++i) {
        int thread = locations[i][0];
        int offset = locations[i][1];
        int size = locations[i][2];
        rw_idx_info_size[i] = size;
        rw_idx.insert(rw_idx.end(),rw_idx_private[thread].begin()+offset,rw_idx_private[thread].begin()+offset+size);
    }
    for (int k = 1; k < graph.n; ++k) {
        rw_idx_info_offset[k] = rw_idx_info_offset[k-1]+rw_idx_info_size[k-1];
    }

//    for (int i = 0; i < graph.n; ++i) {
//        cout<<rw_idx_info_size[i]<<" ";
//    }
//    cout<<endl;
//    for (int i = 0; i < graph.n; ++i) {
//        cout<<rw_idx_info_offset[i]<<" ";
//    }
//    cout<<endl;
//
//    for (int i = 0; i < graph.n; ++i) {
//        cout<<"node "<<i<<endl;
//        for (int j = rw_idx_info_offset[i]; j < rw_idx_info_offset[i]+rw_idx_info_size[i]; ++j) {
//            cout<<rw_idx[j]<<" ";
//        }
//        cout<<endl;
//    }

//    for (auto each:rw_idx_info)
//        cout << each.first<<" "<< each.second<<endl;


    {
        serialize_idx(rwpath);
    }
//#ifdef linux
//    cerr << "Memory usage (MB):" << getMemory()<< endl << endl;
//#endif
}

void backward_push(int target, double rmax,double init_residual=1){
    bwd_idx.first.clean();
    bwd_idx.second.clean();
    vector<int> q;
    q.reserve(graph.n);
    q.push_back(-1);
    unsigned long left = 1;
    q.push_back(target);
    bwd_idx.second.insert(target, init_residual);
    double myeps = rmax;//config.rmax;
    while (left < (int) q.size()) {
        int v = q[left];
        left++;
        double v_residue = bwd_idx.second[v];
        bwd_idx.second[v] = 0;
        if(!bwd_idx.first.exist(v))
            bwd_idx.first.insert( v, v_residue * graph.alpha);
        else
            bwd_idx.first[v] += v_residue * graph.alpha;

        double residual = ((1.0 - graph.alpha) * v_residue);


        for (int next : graph.g[v]) {
            // total_push++;
            double oldep;
            int cnt = graph.g[next].size();
            if( !bwd_idx.second.exist(next) ){
                oldep = 0;
                bwd_idx.second.insert( next,  residual/cnt);
            }
            else{
                oldep = bwd_idx.second[next];
                bwd_idx.second[next] += residual/cnt;
            }
            //if a node's' current residual is small, but next time it got a large residual, it can still be added into forward list
            //so this is correct
            if (oldep <= myeps && bwd_idx.second[next] > myeps) {
                q.push_back(next);
            }
        }
    }

}

int get_dmax(int target, int cluster){
    int dmax=0;
    for(int leaf:level1cluster[cluster]){
        if (leaf==target)
            continue;
        dmax = graph.deg[leaf]>dmax ? graph.deg[leaf]:dmax;
    }
    return dmax;
}


void build_dnpr(int iters = 20){
    double start = omp_get_wtime();
    pr.resize(graph.n,0);
    vector<float> residuals(graph.n, 0);

    for (int i = 0; i < graph.n; ++i) {
        residuals[i] = graph.deg[i]/2.0/graph.m;
    }
    vector<float> new_residuals(graph.n, 0);
    double r_sum = 1;
    uint32_t num_iter = 0;
    while (num_iter<iters) {
        for (int id = 0; id < graph.n; ++id) {
            int degree = graph.deg[id];
            double alpha_residual = graph.alpha * residuals[id];
            pr[id] += alpha_residual;
            r_sum -= alpha_residual;
            double increment = (residuals[id] - alpha_residual) / degree;
            residuals[id] = 0;
            for (int nid:graph.g[id]) {
                new_residuals[nid] += increment;
            }
        }
        residuals.swap(new_residuals);
        num_iter ++ ;
    }
    cout<<omp_get_wtime()-start<<endl;
    std::ofstream ofs(prpath);
    boost::archive::binary_oarchive oa(ofs);
    oa << pr;
}
void build_bwdpush(){
    vector<int> leaf2cluster(graph.n,0);
    for (int i = 0; i < level1cluster.size(); ++i) {
        for(auto leaf:level1cluster[i])
            leaf2cluster[leaf]=i;
    }
    // setting tau

    double tau = graph.tau;
//    double tau = pow(2*graph.k*log(graph.n)*2*graph.m,-1/3.0);
    // collect the node with dnpr>tau
    vector<int> targets;
    unsigned long long size=0;
    targets.reserve(graph.n);
    bwd_idx_info_offset.reserve(graph.n);
    for (int i = 0; i < graph.n; ++i) {
        if (pr[i]>tau){
            targets.push_back(i);
            size+=level1cluster[leaf2cluster[i]].size();
//            cout<<size<<endl;
            bwd_idx_info_offset.push_back(level1cluster[leaf2cluster[i]].size());
        }

    }
//    cout<<"stored target: "<<targets.size()<<endl;
    if (targets.empty()){
        cout<<"no target"<<endl;
        return;
    }
    double start = omp_get_wtime();
    bwd_idx_in_cluster.reserve(size);
    // set error
    double maxval= 1 - pow(graph.n, -2 * graph.epR);
    double minval= 1 - exp(-2 * graph.epR);
    double deltap = graph.n*tau/graph.k;
//    double deltap = tau*graph.dbar;
//    cout<<"scaled delta: "<<deltap<<" 1/10/k: "<<graph.dbar * graph.delta<<endl;
    double epsilonp = min(maxval, max(1 - pow(2 * deltap / exp(1), graph.epR), minval));
    double err = epsilonp*deltap;
    // for each target perform bwdpush
    bwd_idx.first.initialize(graph.n);
    bwd_idx.second.initialize(graph.n);
    double total_time = 0;
    for (int t : targets) {
        int cluster =leaf2cluster[t];
        int dmax = get_dmax(t,cluster);
        double rmax = err/dmax;
//        double start = omp_get_wtime();
        backward_push(t,rmax);
        for(auto source:level1cluster[cluster]){
            if(!bwd_idx.first.exist(source))
                bwd_idx_in_cluster.push_back(0);
            else
                bwd_idx_in_cluster.push_back(graph.deg[source]*bwd_idx.first[source]);
        }
    }
    cout<<omp_get_wtime()-start<<endl;
    string file_name = rwpath+".bwdidx";
    std::ofstream ofs(file_name);
    boost::archive::binary_oarchive oa(ofs);
    oa << bwd_idx_in_cluster;

    string info_name = rwpath+".bwdidx.offset";
    std::ofstream info_ofs(info_name);
    boost::archive::binary_oarchive info_oa(info_ofs);
    info_oa << bwd_idx_info_offset;
}


bool exists_test(const std::string &name) {
    ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    }
    else {
        f.close();
        return false;
    }
}

void assert_file_exist(string desc, string name) {

    if (!exists_test(name)) {
        cerr << desc << " " << name << " not find " << endl;
        exit(1);
    }
}

inline void deserialize_idx(){
    string file_name = rwpath+".idx";
    assert_file_exist("index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    ia >> rw_idx;

    string info_name = rwpath+".info.offset";
    assert_file_exist("index file", file_name);
    std::ifstream info_ifs(info_name);
    boost::archive::binary_iarchive info_ia(info_ifs);
    info_ia >> rw_idx_info_offset;

    string size_name = rwpath+".info.size";
    assert_file_exist("index file", size_name);
    std::ifstream size_ifs(size_name);
    boost::archive::binary_iarchive size_ia(size_ifs);
    size_ia >> rw_idx_info_size;

}

inline void deserialize_pr(){
    string file_name = prpath;
    assert_file_exist("dnpr index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    ia >> pr;
}

inline void deserialize_bwd(){
    double tau = graph.tau;
    // collect the node with dnpr>tau
    vector<int> targets;
    targets.reserve(graph.n);
    for (int i = 0; i < graph.n; ++i) {
        if (pr[i]>tau)
            targets.push_back(i);
    }
    if (targets.empty()){
        return;
    }
    string file_name = rwpath+".bwdidx";
    assert_file_exist("index file", file_name);
    std::ifstream ifs(file_name);
    boost::archive::binary_iarchive ia(ifs);
    ia >> bwd_idx_in_cluster;

    string info_name = rwpath+".bwdidx.offset";
    assert_file_exist("index file", file_name);
    std::ifstream info_ifs(info_name);
    boost::archive::binary_iarchive info_ia(info_ifs);
    info_ia >> bwd_idx_info_offset;
    int offset=0;
    for (int i = 0; i < targets.size(); ++i) {
        int t = targets[i];
        int size = bwd_idx_info_offset[i];
        bwd_idx_info.insert({t,make_pair(offset,size)});
        offset += size;
    }
}

void init_container(){
    fwd_idx.resize(thread_num);
    virtual_fwd_idx.resize(thread_num);
    queues.resize(thread_num);
    for (int i = 0; i < thread_num; ++i) {
        fwd_idx[i].first.initialize(graph.n);
        fwd_idx[i].second.initialize(graph.n);
        virtual_fwd_idx[i].first.initialize(graph.n);
        virtual_fwd_idx[i].second.initialize(graph.n);
        queues[i].reserve(graph.n);
    }

//    rsums.resize(graph.n);
    // initialize stack info
//    fwd_stack.resize(graph.max_level);
    leaf2id_stack.resize(graph.max_level);
    fwd_stack_top=0;
}

void forward_push(vector<int> &q, unsigned long left,
                  double &rsum, double rmax,int tid) {
    Fwdidx &fwdidx = fwd_idx[tid];
    double myeps = rmax;//config.rmax;
    while (left < (int) q.size()) {
        int v = q[left];
        left++;
        double v_residue = fwdidx.second[v];
        fwdidx.second[v] = 0;
        if(!fwdidx.first.exist(v))
            fwdidx.first.insert( v, v_residue * graph.alpha);
        else
            fwdidx.first[v] += v_residue * graph.alpha;

        int out_neighbor = graph.deg[v];
        rsum -=v_residue*graph.alpha;
//        assert(out_neighbor > 0);
        double avg_push_residual = ((1.0 - graph.alpha) * v_residue) / out_neighbor;

        for (int next : graph.g[v]) {
            // total_push++;
            double oldep;
            if( !fwdidx.second.exist(next) ){
                oldep = 0;
                fwdidx.second.insert( next,  avg_push_residual);
            }
            else{
                oldep = fwdidx.second[next];
                fwdidx.second[next] += avg_push_residual;
            }
            //if a node's' current residual is small, but next time it got a large residual, it can still be added into forward list
            //so this is correct
            if (oldep/graph.deg[next] <= myeps && fwdidx.second[next]/graph.deg[next] > myeps) {
                q.push_back(next);
            }
        }
    }
    if(verbose){
        double correct = 0;
        for(long i=0; i < fwdidx.second.occur.m_num; i++){
            int source = fwdidx.second.occur[i];
            correct += fwdidx.second[source];
        }
        assert(abs(rsum-correct)<10e-6);
    }
}



void update_cache_virtual_fwd(int leaf, int tid){
    Fwdidx &vfwdidx = virtual_fwd_idx[tid];
    Fwdidx &fwdidx = fwd_idx[tid];
//    if(visual_mode == INTERACT_MODE){
    int pos = (*leaf2id)[leaf];

    for(long i=0; i < fwdidx.first.occur.m_num; i++){
        int source = fwdidx.first.occur[i];
        double reserve = fwdidx.first[source];
        // insert for reserve_cache
//            block.pdst[i] = source;
//            block.p[i] = reserve;
        if (reserve==0)
            continue;
        // update virtual_fwd_idx
        if( !vfwdidx.first.exist(source) )
            vfwdidx.first.insert(source, reserve);
        else
            vfwdidx.first[source] += reserve;
//        if (verbose)
//            cout<<"virtual_fwd_idx p"<<source<<": "<<virtual_fwd_idx.first[source]<<endl;
    }
//        block.rdst.resize(fwdidx.second.occur.m_num);
//        block.r.resize(fwdidx.second.occur.m_num);
    for(long i=0; i < fwdidx.second.occur.m_num; i++){
        int source = fwdidx.second.occur[i];
        double residual = fwdidx.second[source];
        // insert for reside_cache
//            block.rdst[i] = source;
//            block.r[i] = residual;
        if (residual==0)
            continue;
        // update virtual_fwd_idx
        if( !vfwdidx.second.exist(source) )
            vfwdidx.second.insert(source, residual);
        else
            vfwdidx.second[source] += residual;
//        if (verbose)
//            cout<<"virtual_fwd_idx r"<<source<<": "<<virtual_fwd_idx.second[source]<<endl;
    }
}




void forward_local_update_linear(int s, const string &Father, double& rsum, double rmax, int tid){
    double init_residual;
    vector<int> *ids;
    Fwdidx &fwdidx = fwd_idx[tid];

    init_residual = graph.deg[s];

    // load reserve and residue info from stack
    // initialize queue and fwd_idx by rescaling the cache results
    // set rmax for leaf forward push
    double myeps = rmax;
    // set rmax for leaf forward push
    double r = init_residual;
    if(r / graph.deg[s] < myeps){
        rsum=r;
        return;
    }

    vector<int> &q = queues[tid];
    q.push_back(-1);
    unsigned long left = 1;
    // make sure fwd_idx is clean before usage
    assert(fwdidx.second.occur.m_num==0);
    assert(fwdidx.first.occur.m_num==0);

//    if (fwd_stack_top==0 || visual_mode==FULL_MODE){
    // initial queue
    q.push_back(s);

    fwdidx.second.insert(s, r);

    // start forward push
    forward_push(q, left, r, myeps,tid);
    // collect the rsum from each leaf source
    rsum=r;
    q.clear();
    q.reserve(graph.n);
//    cout<<"nnz: "<<fwdidx.second.occur.m_num<<endl;
//    fwd_idx.first.clean();
//    fwd_idx.second.clean();
}

void forward_local_update_all_leaf(const string &S, const string &Father,
                                   double& rsum, double rmax, int tid){
    vector<int> *ids;
//    if(fwd_stack_top==0)//    else{
//        ids = &leaf2id_stack[fwd_stack_top-1];
//        init_residual = (double)super2leaf[Father].size()/(double)super2leaf[S].size();
//    }
    // load reserve and residue info from stack
//    double nnz = 0;
    for(int leaf:super2leaf[S]){
        double myeps = sqrt(graph.deg[leaf])*rmax; // d(leaf)/dsum fraction of rmax;
//        /graph.g[next].size() > myeps
        // set rmax for leaf forward push
        double r = graph.deg[leaf];
        vector<int> &q = queues[tid];
        q.push_back(-1);
        unsigned long left = 1;
        // make sure fwd_idx is clean before usage
        pair<iMap<double>, iMap<double>> &fwdidx = fwd_idx[tid];
        assert(fwdidx.second.occur.m_num == 0);
        assert(fwdidx.first.occur.m_num == 0);
        // set rmax for leaf forward push

//        if(fwd_stack_top==0){
        // initial queue
        q.push_back(leaf);
        // initial rsum
        // rsums[leaf] = r;
        // initial fwd_idx
        fwdidx.second.insert(leaf, r);
        // start forward push
        forward_push(q, left, r, myeps,tid);
        // collect the rsum from each leaf source
        rsum+=r;
        // store fwd_idx info into cache and virtual_fwd_idx
        update_cache_virtual_fwd(leaf,tid);
//        total_rsum += rsum;

        // clean fwd_idx
        q.clear();
        q.reserve(graph.n);
        fwdidx.first.clean();
        fwdidx.second.clean();
    }
//    nnz += virtual_fwd_idx[tid].second.occur.m_num;
//    cout<<"nnz r: "<<nnz<<" leaf size: "<<super2leaf[S].size()<<" ratio: "<<nnz/super2leaf[S].size()<<endl;
    // check rsum is correct here
//    if(verbose)
//        assert(is_checked_rsum(rsum,tid));

}

void forward_local_update_all_virtual(const string &S, const string &Father,
                                      double& rsum, double myeps, int dsum, int tid){
    double init_residual;
    vector<int> *ids;
    double leafsum = super2leaf[S].size();
    rsum = dsum/leafsum;
    vector<int> &q = queues[tid];
    q.push_back(-1);
    unsigned long left = 1;
    // make sure fwd_idx is clean before usage
    Fwdidx &vfwdidx = virtual_fwd_idx[tid];
    assert(vfwdidx.second.occur.m_num == 0);
    assert(vfwdidx.first.occur.m_num == 0);
//    double nnz = 0;

    if (super2leaf[S].size()>1/myeps){
        for(int leaf:super2leaf[S]){
            // initialize queue and fwd_idx by rescaling the cache results
            vfwdidx.second.insert(leaf, graph.deg[leaf]/leafsum);
        }
    } else{
        for(int leaf:super2leaf[S]){
            // initialize queue and fwd_idx by rescaling the cache results
            q.push_back(leaf);
            vfwdidx.second.insert(leaf, graph.deg[leaf]/leafsum);
        }
        // start forward push
        while (left < (int) q.size()) {
            int v = q[left];
            left++;
            double v_residue = vfwdidx.second[v];
            vfwdidx.second[v] = 0;
            if(!vfwdidx.first.exist(v))
                vfwdidx.first.insert( v, v_residue * graph.alpha);
            else
                vfwdidx.first[v] += v_residue * graph.alpha;

            int out_neighbor = graph.deg[v];
            rsum -=v_residue*graph.alpha;
//            assert(out_neighbor > 0);
            double avg_push_residual = ((1.0 - graph.alpha) * v_residue) / out_neighbor;

            for (int next : graph.g[v]) {
                // total_push++;
                double olderes;
                if( !vfwdidx.second.exist(next) ){
                    olderes = 0;
                    vfwdidx.second.insert( next,  avg_push_residual);
                }
                else{
                    olderes = vfwdidx.second[next];
                    vfwdidx.second[next] += avg_push_residual;
                }
                //if a node's' current residual is small, but next time it got a large residual, it can still be added into forward list
                //so this is correct
                if (olderes / graph.deg[next] <= myeps && vfwdidx.second[next] / graph.deg[next] > myeps) {
                    q.push_back(next);
                }
            }
        }
    }



//    total_rsum += rsum;
//    nnz += vfwdidx.second.occur.m_num;
    // clean fwd_idx
    q.clear();
    q.reserve(graph.n);
}


void compute_ppr_with_fwdidx(int S, const string &Father,unsigned long long &num_random_walk,
                             double omega, double check_rsum, int tid, vector<bool> &insideLeaf,
                             vector<int> &spprs,vector<int> &spprt,vector<double> &spprv,
                             unsigned int &seed){
    // INFO("rsum is:", check_rsum);
    Fwdidx &fwdidx = fwd_idx[tid];
    unordered_map<int,double> pprdegs;
    for(long i=0; i < fwdidx.first.occur.m_num; i++){
        int leaf = fwdidx.first.occur[i];
        double val = fwdidx.first[leaf];
        if (insideLeaf[leaf])
            pprdegs[leaf] += val;
    }
//    if(omega * check_rsum < 1){
//        for(auto pair:pprdegs){
//            int T = pair.first;
//            spprs.push_back(S);
//            spprt.push_back(T);
//            spprv.push_back(pair.second);
//        }
//        return;
//    }

    // INFO(num_random_walk);
    double ideal_walk_number = omega * check_rsum;
    double incre = check_rsum/ideal_walk_number;
    //rand walk online
    if (graph.isOnline){
        fwdidx.second.occur.Sort();
        for(long i=0; i < fwdidx.second.occur.m_num; i++){
//            if (num_random_walk>=ideal_walk_number)
//                break;
            int source = fwdidx.second.occur[i];
            double residual = fwdidx.second[source];

//            if (residual<threshold)
            if (residual==0)
                continue;
//            nnz++;
            unsigned long num_s_rw = ceil(residual*omega);
//            num_total_rw += num_s_rw;
//            num_rw_level += num_s_rw;
            num_random_walk+=num_s_rw;


            for(unsigned long j=0; j<num_s_rw; j++){
                int des = random_walk(source, seed);
                if(insideLeaf[des]) pprdegs[des] += incre;
            }
        }
    }else{
        fwdidx.second.occur.Sort();
//        unsigned long olcnt = 0;
        for(long i=0; i < fwdidx.second.occur.m_num; i++){
//            if (num_random_walk>=ideal_walk_number)
//                break;
            int source = fwdidx.second.occur[i];
            double residual = fwdidx.second[source];

//            if(residual<threshold)
            if (residual==0)
                continue;
//            usennz++;
            unsigned long num_s_rw = ceil(residual*omega);
//            num_total_rw += num_s_rw;
//            num_rw_level += num_s_rw;
            num_random_walk += num_s_rw;
            used_rw_size += num_s_rw;
            idx_rw_size += rw_idx_info_size[source];
            //for each source node, get rand walk destinations from previously generated idx or online rand walks
            if(num_s_rw > rw_idx_info_size[source]){ //if we need more destinations than that in idx, rand walk online
//                olcnt+=num_s_rw-rw_idx_info_size[source];
                for(unsigned long k=0; k<rw_idx_info_size[source]; k++){
                    int des = rw_idx[rw_idx_info_offset[source] + k];
                    if(insideLeaf[des]) pprdegs[des] += incre;
                }
//                num_hit_idx += rw_idx_info_size[source];

                for(unsigned long j=0; j < num_s_rw-rw_idx_info_size[source]; j++){ //rand walk online
                    int des = random_walk(source,seed);
                    if(insideLeaf[des]) pprdegs[des] += incre;
                }
            }else{ // using previously generated idx is enough
                for(unsigned long k=0; k<num_s_rw; k++){
                    int des = rw_idx[rw_idx_info_offset[source] + k];
                    if(insideLeaf[des]) pprdegs[des] += incre;
                }
//                num_hit_idx += num_s_rw;
            }
        }
//        cout << "ol rw: "<<olcnt<<endl;
    }
    for(auto pair:pprdegs){
        int T = pair.first;
        spprs.push_back(S);
        spprt.push_back(T);
        spprv.push_back(pair.second);
    }

}

void compute_super_ppr_with_fwdidx(const string &S, const string &Father,
                                   unsigned long long &num_random_walk, double omega,
                                   double check_rsum, int tid,vector<int> &leaf2super,
                                   vector<int> &spprs,vector<int> &spprt,vector<double> &spprv,
                                   unsigned int &seed){
    // INFO("rsum is:", check_rsum);
    int compo,level;
    int Sid = extract_supernode(S,compo,level);
    double ideal_walk_number = omega * check_rsum;
    Fwdidx &vfwdidx = virtual_fwd_idx[tid];
    unordered_map<int,double> pprdegs;
//    double start = omp_get_wtime();
    for(long i=0; i < vfwdidx.first.occur.m_num; i++){
        int leaf = vfwdidx.first.occur[i];
        int Tid = leaf2super[leaf];
        if (Tid==-1) continue;
        double val = vfwdidx.first[leaf];
        pprdegs[Tid] += val;
    }
//    fwdTime[2] += omp_get_wtime()-start;

//    if( ideal_walk_number <= 1){
//        for(auto pair:pprdegs){
//            int Tid = pair.first;
//            string Tname = "c"+to_string(compo)+"_l"+to_string(level)+"_"+to_string(Tid);
//            spprs.push_back(Sid);
//            spprt.push_back(Tid);
//            spprv.push_back(pair.second / (double)super2leaf[Tname].size());
//        }
//        return;
//    }
//    start = omp_get_wtime();
    // INFO(num_random_walk);
    double incre = check_rsum/(double)ideal_walk_number;
    //rand walk online
    if (graph.isOnline){
        vfwdidx.second.occur.Sort();
        for(long i=0; i < vfwdidx.second.occur.m_num; i++){
//            if (num_random_walk>=ideal_walk_number)
//                break;
            int source = vfwdidx.second.occur[i];
            double residual = vfwdidx.second[source];
//            if (residual<threshold)
            if (residual==0)
                continue;
//            nnz++;
            unsigned long num_s_rw = ceil(residual*omega);
//            num_total_rw += num_s_rw;
//            num_rw_level += num_s_rw;
            num_random_walk+=num_s_rw;


            for(unsigned long j=0; j<num_s_rw; j++){
                int des = random_walk(source,seed);
                int Tid = leaf2super[des];
                if (Tid==-1) continue;
                pprdegs[Tid] += incre;
            }
        }
    }else{
        vfwdidx.second.occur.Sort();
        for(long i=0; i < vfwdidx.second.occur.m_num; i++){
//            if (num_random_walk>=ideal_walk_number)
//                break;
            int source = vfwdidx.second.occur[i];
            double residual = vfwdidx.second[source];
//            if(residual<threshold)
            if (residual==0)
                continue;
//            usennz ++;
            unsigned long num_s_rw = ceil(residual*omega);
//            num_total_rw += num_s_rw;
//            num_rw_level += num_s_rw;
            num_random_walk += num_s_rw;
            used_rw_size += num_s_rw;
            idx_rw_size += rw_idx_info_size[source];
            //for each source node, get rand walk destinations from previously generated idx or online rand walks
            if(num_s_rw > rw_idx_info_size[source]){ //if we need more destinations than that in idx, rand walk online
                for(unsigned long k=0; k<rw_idx_info_size[source]; k++){
                    int des = rw_idx[rw_idx_info_offset[source] + k];
                    int Tid = leaf2super[des];
                    if (Tid==-1) continue;
                    pprdegs[Tid] += incre;
                }
//                num_hit_idx += rw_idx_info_size[source];

                for(unsigned long j=0; j < num_s_rw-rw_idx_info_size[source]; j++){ //rand walk online
                    int des = random_walk(source,seed);
                    int Tid = leaf2super[des];
                    if (Tid==-1) continue;
                    pprdegs[Tid] += incre;
                }
            }else{ // using previously generated idx is enough
                for(unsigned long k=0; k<num_s_rw; k++){
                    int des = rw_idx[rw_idx_info_offset[source] + k];
                    int Tid = leaf2super[des];
                    if (Tid==-1) continue;
                    pprdegs[Tid] += incre;
                }
//                num_hit_idx += num_s_rw;
            }
        }
    }
//    fwdTime[3] += omp_get_wtime()-start;

    for(auto pair:pprdegs){
        int Tid = pair.first;
        string Tname = "c"+to_string(compo)+"_l"+to_string(level)+"_"+to_string(Tid);
        spprs.push_back(Sid);
        spprt.push_back(Tid);
        spprv.push_back(pair.second / (double)super2leaf[Tname].size());
    }

}


void foratp_super_pprdeg(const string &S, const string &Father,
                         double rmax, double omega, int dsum, int level, int tid, vector<int> &leaf2super,
                         vector<int> &spprs, vector<int> &spprt, vector<double> &spprv,
                         unsigned int &seed){
//    double start = omp_get_wtime();
    double rsum=0;
    //forward propagation, obtain reserve and residual
    forward_local_update_all_virtual(S, Father, rsum, rmax, dsum, tid);

//    fwdTime[rwTime.size()-level] += omp_get_wtime()-start;
    //reserve mem for superPPRDeg
    Fwdidx & vfwdidx = virtual_fwd_idx[tid];
    spprs.reserve(vfwdidx.first.occur.m_num);
    spprt.reserve(vfwdidx.first.occur.m_num);
    spprv.reserve(vfwdidx.first.occur.m_num);
//    cout<<"super fwd time:"<<omp_get_wtime()-start<<endl;
    if (verbose){
//        cout << "fwd time:"<<timeBy(start)<<endl;
        cout <<"rmax: "<<rmax <<endl;
    }
//    graph.rmaxs += rmax;
//    graph.rsums += rsum;
//    start = omp_get_wtime();
    unsigned long long num_random_walk=0;
    compute_super_ppr_with_fwdidx(S,Father,num_random_walk,omega, rsum,tid,leaf2super,spprs,spprt,spprv,seed);
//    cout << "omega * check_rsum: "<<omega * rsum<<endl;
//    cout<<"super rw time:"<<omp_get_wtime()-start<<endl;
//    cout<<"walk number: "<<num_random_walk<<endl;
//    rwTime[rwTime.size()-level] += omp_get_wtime()-start;
    // clear virtual_fwd_idx
    vfwdidx.first.clean();
    vfwdidx.second.clean();
    if (verbose){
        cout << "rsum: "<<rsum<<endl;
        cout << "omega: "<<omega<<endl;
        cout << "rsum*omega: "<<rsum*omega<<endl;
//        cout<<"walk time: "<<timeBy(start)<<endl;
//        cout << "nnz: "<<nnz<<endl;
        cout<<"walk number: "<<num_random_walk<<endl;
        cout<<"--------------"<<endl;
    }
//    nnz=0;
}


void fora_super_pprdeg(const string &S, const string &Father,
                       double rmax, double omega, int level, int tid, vector<int> &leaf2super,
                       vector<int> &spprs, vector<int> &spprt, vector<double> &spprv,
                       unsigned int &seed){
//    double start = omp_get_wtime();
    double rsum=0;
    //forward propagation, obtain reserve and residual
    forward_local_update_all_leaf(S, Father, rsum, rmax,  tid);

//    fwdTime[rwTime.size()-level] += omp_get_wtime()-start;
    //reserve mem for superPPRDeg
    Fwdidx & vfwdidx = virtual_fwd_idx[tid];
    spprs.reserve(vfwdidx.first.occur.m_num);
    spprt.reserve(vfwdidx.first.occur.m_num);
    spprv.reserve(vfwdidx.first.occur.m_num);
//    cout<<"super fwd time:"<<omp_get_wtime()-start<<endl;
    if (verbose){
//        cout << "fwd time:"<<timeBy(start)<<endl;
        cout <<"rmax: "<<rmax <<endl;
    }
//    graph.rmaxs += rmax;
//    graph.rsums += rsum;
//    start = omp_get_wtime();
    unsigned long long num_random_walk=0;
    compute_super_ppr_with_fwdidx(S,Father,num_random_walk,omega, rsum,tid,leaf2super,spprs,spprt,spprv,seed);
    int size = super2leaf[S].size();
    for (double & i : spprv) {
        i/=(double)size;
    }
//    rwTime[rwTime.size()-level] += omp_get_wtime()-start;
    // clear virtual_fwd_idx
    vfwdidx.first.clean();
    vfwdidx.second.clean();

}

void taupush_super_pprdeg(const string &S, const string &Father,
                          double rmax, int dsum, int level, int tid, vector<int> &leaf2super,
                          vector<int> &spprs, vector<int> &spprt, vector<double> &spprv,
                          unsigned int &seed){
//    double start = omp_get_wtime();
    double rsum=0;
    //forward propagation, obtain reserve and residual
    forward_local_update_all_virtual(S, Father, rsum, rmax, dsum, tid);


//    fwdTime[rwTime.size()-level] += omp_get_wtime()-start;
    //reserve mem for superPPRDeg
    Fwdidx & vfwdidx = virtual_fwd_idx[tid];



    spprs.reserve(vfwdidx.first.occur.m_num);
    spprt.reserve(vfwdidx.first.occur.m_num);
    spprv.reserve(vfwdidx.first.occur.m_num);
//    cout<<"super fwd time:"<<omp_get_wtime()-start<<endl;
    if (verbose){
//        cout << "fwd time:"<<timeBy(start)<<endl;
        cout <<"rmax: "<<rmax <<endl;
    }
    int compo;
    int Sid = extract_supernode(S,compo,level);
    unordered_map<int,double> pprdegs;
    for(long i=0; i < vfwdidx.first.occur.m_num; i++){
        int leaf = vfwdidx.first.occur[i];
        int Tid = leaf2super[leaf];
        if (Tid==-1) continue;
        double val = vfwdidx.first[leaf];
        pprdegs[Tid] += val;
    }
    for(auto pair:pprdegs){
        int Tid = pair.first;
        string Tname = "c"+to_string(compo)+"_l"+to_string(level)+"_"+to_string(Tid);
        spprs.push_back(Sid);
        spprt.push_back(Tid);
        spprv.push_back(pair.second / (double)super2leaf[Tname].size());
    }
    // clear virtual_fwd_idx
    vfwdidx.first.clean();
    vfwdidx.second.clean();

//    nnz=0;
}

void power_iteration(vector<double> &pi, vector<double> &residuals, double &r_sum, double l1_error=1e-9){
    vector<double> new_residuals(graph.n, 0);
    for (uint32_t num_iter = 0; r_sum > l1_error; ++num_iter) {
        for (int id = 0; id < graph.n; ++id) {
            int degree = graph.deg[id];
            double alpha_residual = graph.alpha * residuals[id];
            pi[id] += alpha_residual;
            r_sum -= alpha_residual;
            double increment = (residuals[id] - alpha_residual) / degree;
            residuals[id] = 0;
            for (int nid:graph.g[id]) {
                new_residuals[nid] += increment;
            }
        }
        residuals.swap(new_residuals);
    }
}

void powiter_super_pprdeg(const string &S, const string &Father,
                          int level, int tid, vector<int> &leaf2super,
                          vector<int> &spprs, vector<int> &spprt, vector<double> &spprv,
                          unsigned int &seed){
//    double start = omp_get_wtime();
    double rsum=0;
    vector<double> pi(graph.n, 0);
    vector<double> residuals(graph.n, 0);
    unordered_map<int,double> pprdegs;
    for(int leaf:super2leaf[S]){
        fill(pi.begin(),pi.end(),0);
        fill(residuals.begin(),residuals.end(),0);
        rsum=graph.deg[leaf];
        residuals[leaf] = graph.deg[leaf];
        power_iteration(pi, residuals, rsum);
        for(long i=0; i < residuals.size(); i++){
            int Tid = leaf2super[i];
            if (Tid==-1) continue;
            double val = pi[i];
            pprdegs[Tid] += val;
        }
    }

//    fwdTime[rwTime.size()-level] += omp_get_wtime()-start;
    int compo;
    int Sid = extract_supernode(S,compo,level);
    int source_size = super2leaf[S].size();
    int size = pprdegs.size();
    spprs.reserve(size);
    spprt.reserve(size);
    spprv.reserve(size);
    for(auto pair:pprdegs){
        int Tid = pair.first;
        string Tname = "c"+to_string(compo)+"_l"+to_string(level)+"_"+to_string(Tid);
        spprs.push_back(Sid);
        spprt.push_back(Tid);
        spprv.push_back(pair.second / (double)super2leaf[Tname].size()/(double)source_size);
    }
}

int static get_indice(int i,int j, int n){
    return i*n+j;
}

void config_superpprdeg(int level, vector<string> &childrens, double &epsilonR, double &delta, int &gamma, vector<int> &leafsize,
                        vector<int> &dsum, double &dmax, vector<int> &leaf2super, double &Delta, double &maxval, double &minval, int &nrow) {
    epsilonR= graph.epR;
    delta= graph.delta;
    gamma= INT_MAX;
    Delta= 0;
    maxval= 1 - pow(graph.n, -2 * epsilonR);
    minval= 1 - exp(-2 * epsilonR);
    nrow= childrens.size();
//    double delta = pow(25,level-1)*graph.delta;
// tau = min_{C\in V_i}|Leaf(C)|
// record |Leaf(C)| for each C
// record d_{sum}(C) for each C
    int id = 0;
    vector<int> &ids = *leaf2id;
    ids.resize(graph.n,-1);
#pragma omp parallel num_threads(thread_num)
//#pragma omp single
#pragma omp for schedule(static) nowait
    for (int i = 0; i < childrens.size(); ++i) {
        string &child = childrens[i];
        int size = super2leaf[child].size();
        leafsize[i] = size;
        if (size < gamma)
            gamma = size;
        vector<int> &leafs = super2leaf[child];
        for(int j=0;j<leafs.size();++j){
            int leaf = leafs[j];
            int _;
            leaf2super[leaf] = extract_supernode(child,_,_);
            int deg = graph.deg[leaf];
            dsum[i] += deg;
            ids[leaf] = id;
            id++;
        }
    }

//    if (isVNFP==1){
//        for (int i = 0; i < childrens.size(); ++i) {
//            Delta += (double)dsum[i]/leafsize[i];
//        }
//    }
//    else {
//        for (int i = 0; i < childrens.size(); ++i) {
//            Delta += (double)dsum[i];
//        }
//    }
    for (int i = 0; i < childrens.size(); ++i) {
        double davg = (double)dsum[i]/leafsize[i];
        Delta += davg;
        dmax = max(dmax,davg);
    }

}

void config_pr_superpprdeg(int level, vector<string> &childrens, double &epsilonR, double &tau, vector<int> &leafsize,
                           vector<int> &dsum, vector<int> &leaf2super, double &maxval, double &minval, int &nrow,vector<string> &large_target) {
    epsilonR= graph.epR;
    double max_tau = tau;
    tau= 0;
    vector<double> local_tau(childrens.size(),0);
    maxval= 1 - pow(graph.n, -2 * epsilonR);
    minval= 1 - exp(-2 * epsilonR);
    nrow= childrens.size();
    int id = 0;
    vector<int> &ids = *leaf2id;
    ids.resize(graph.n,-1);
#pragma omp parallel num_threads(thread_num)
//#pragma omp single
#pragma omp for schedule(static) nowait
    for (int i = 0; i < childrens.size(); ++i) {
        string &child = childrens[i];
        int size = super2leaf[child].size();
        leafsize[i] = size;
        vector<int> &leafs = super2leaf[child];
        double& supertau = local_tau[i];
        for(int j=0;j<leafs.size();++j){
            int leaf = leafs[j];
            supertau += pr[leaf];
            int _;
            leaf2super[leaf] = extract_supernode(child,_,_);
            dsum[i] += graph.deg[leaf];
            ids[leaf] = id;
            id++;
        }
        supertau/=leafs.size();
    }
    if (isBPSN){
        priority_queue<pair<float,string>> sort_tau;
        for (int i = 0; i < childrens.size(); ++i) {
            string &child = childrens[i];
            sort_tau.push(pair<float, string>(local_tau[i], child));
        }
        while (sort_tau.top().first>max_tau){
            large_target.push_back(sort_tau.top().second);
            sort_tau.pop();
        }
        tau = sort_tau.top().first;
    }
    else{
        for(double each:local_tau){
            tau = max(each,tau);
        }
    }

}

void allpair_super_pprdeg(const string &supernode, int level,
                          vector<string> &childrens, vector<int> &node2id){

    double start,epsilonR,delta,Delta,epsilonp, maxval, minval;
    double dmax=0;
    int gamma,nrow;
    vector<int> leafsize(childrens.size());
    vector<int> dsum(childrens.size());
    vector<int> leaf2super(graph.n,-1);
    config_superpprdeg(level, childrens, epsilonR, delta, gamma, leafsize, dsum, dmax, leaf2super, Delta, maxval, minval, nrow);
    start = omp_get_wtime();
    int childsize = childrens.size();
    double mint=10000,maxt=0;
#pragma omp parallel num_threads(thread_num)
//#pragma omp single
#pragma omp for schedule(static) nowait
    for (int i = 0; i < childrens.size(); ++i) {
        double rmax = 0;
        double omega = 0;
        double s1 = omp_get_wtime();
        assert(omp_get_num_threads()==thread_num);
        string &child = childrens[i];
        vector<int> spprs, spprt;
        vector<double> spprv;
        // parameter settings for each ss superppr
        // set deltap for current source supernode as the average degree of its leaf nodes
        double deltap = graph.dbar * delta;
//        double deltap = (double) dsum[i] / leafsize[i] * delta;
        // set the corresponding epsilonp for current source supernode
        epsilonp = min(maxval, max(1 - pow(2 * deltap / exp(1), epsilonR), minval));
        // rmax and omega settings for FORA
        if (isFORASN){
            rmax = epsilonp * sqrt(gamma * Delta * deltap / (2.0 * graph.m) / (2 + 2 * epsilonp / 3.0) / log(1 / graph.pfail));
            omega = (2 + 2 * epsilonp / 3.0) * log(1 / graph.pfail) / epsilonp / epsilonp / deltap / gamma;
        }
        else{
            rmax = epsilonp * sqrt( deltap / (2.0 * graph.m) / (2 + 2 * epsilonp / 3.0) / log(1 / graph.pfail));
            omega = (2 + 2 * epsilonp / 3.0) * log(1 / graph.pfail) / epsilonp / epsilonp / deltap;
        }

//        }

        if (verbose) {
            cout << "Delta: " << Delta << endl;
            cout << "tau: " << gamma << endl;
            cout << "deltap: " << deltap << endl;
            cout << "epsilonp: " << epsilonp << endl;
            cout << "deltap*epsilonp: " << epsilonp*deltap << endl;
        }
        unsigned int seed = (i+1)*(omp_get_thread_num()+1);
        if (isFORASN)
            foratp_super_pprdeg(child, supernode, rmax, omega, dsum[i], level, omp_get_thread_num(), leaf2super, spprs,
                                spprt, spprv, seed);
        else
            fora_super_pprdeg(child, supernode, rmax, omega, level, omp_get_thread_num(), leaf2super, spprs,spprt, spprv, seed);
        long leng = spprs.size();
        int idx;

        for (int i = 0; i < leng; ++i) {
            idx = get_indice(node2id[spprs[i]], node2id[spprt[i]], nrow);
            PPRDeg[idx] = spprv[i];
        }


    }
    timeElasped += omp_get_wtime()-start;
//    levelTime[levelTime.size()-level] += omp_get_wtime()-start;
}


void allpair_super_pprdeg_powiter(const string &supernode, int level,
                                  vector<string> &childrens, vector<int> &node2id){
    double start;
    int nrow=childrens.size();;
    vector<int> leaf2super(graph.n,-1);
#pragma omp parallel num_threads(thread_num)
//#pragma omp single
#pragma omp for schedule(static) nowait
    for (int i = 0; i < childrens.size(); ++i) {
        string &child = childrens[i];
        int _;
        for(int leaf:super2leaf[child]){
            leaf2super[leaf] = extract_supernode(child,_,_);
        }
    }
    start = omp_get_wtime();
#pragma omp parallel num_threads(thread_num)
//#pragma omp single
#pragma omp for schedule(static) nowait
    for (int i = 0; i < childrens.size(); ++i) {
        double s1 = omp_get_wtime();
        assert(omp_get_num_threads()==thread_num);
        string &child = childrens[i];
        vector<int> spprs, spprt;
        vector<double> spprv;

        unsigned int seed = (i+1)*(omp_get_thread_num()+1);
        powiter_super_pprdeg(child, supernode,level, omp_get_thread_num(), leaf2super, spprs, spprt,
                             spprv, seed);
        long leng = spprs.size();
        int idx;
        for (int i = 0; i < leng; ++i) {
            idx = get_indice(node2id[spprs[i]], node2id[spprt[i]], nrow);
            PPRDeg[idx] = spprv[i];
        }

    }
    timeElasped += omp_get_wtime()-start;
//    levelTime[levelTime.size()-level] += omp_get_wtime()-start;
}

void allpair_super_pprdeg_pr(const string &supernode, int level,
                             vector<string> &childrens, vector<int> &node2id){

    double start,epsilonR,epsilonp, maxval, minval,tau;
    int nrow;
    if (isBPSN){
        tau = graph.tau;
    } else{
        tau = 0;
    }

    vector<int> leafsize(childrens.size());
    vector<int> dsum(childrens.size());
    vector<int> leaf2super(graph.n,-1);
    vector<string> large_target;
    config_pr_superpprdeg(level, childrens, epsilonR, tau, leafsize, dsum, leaf2super, maxval, minval, nrow,large_target);
    start = omp_get_wtime();
    int childsize = childrens.size();
    double mint=10000,maxt=0;
    double total_rinit = 0;
#pragma omp parallel num_threads(thread_num)
//#pragma omp single
#pragma omp for schedule(static) nowait
    for (int i = 0; i < childrens.size(); ++i) {
        double rmax = 0;
        double s1 = omp_get_wtime();
        assert(omp_get_num_threads()==thread_num);
        string &child = childrens[i];
        vector<int> spprs, spprt;
        vector<double> spprv;
        // parameter settings for each ss superppr
        // set deltap for current source supernode as the average degree of its leaf nodes
        double deltap = graph.dbar * graph.delta;
//        double deltap = (double) dsum[i] / leafsize[i] * delta;
        // set the corresponding epsilonp for current source supernode
        epsilonp = min(maxval, max(1 - pow(2 * deltap / exp(1), epsilonR), minval));
        // rmax settings
        rmax = epsilonp*deltap/tau/graph.m/2;


        if (verbose) {
            cout << "tau: " << tau << endl;
            cout << "deltap: " << deltap << endl;
            cout << "epsilonp: " << epsilonp << endl;
            cout << "deltap*epsilonp: " << epsilonp*deltap << endl;
        }
        unsigned int seed = (i+1)*(omp_get_thread_num()+1);
        taupush_super_pprdeg(child, supernode, rmax, dsum[i], level, omp_get_thread_num(), leaf2super, spprs, spprt,
                             spprv, seed);
        long leng = spprs.size();
        int idx;

        for (int i = 0; i < leng; ++i) {
            idx = get_indice(node2id[spprs[i]], node2id[spprt[i]], nrow);
            PPRDeg[idx] = spprv[i];
        }
        assert(large_target.empty());
    }
    timeElasped += omp_get_wtime()-start;
//    levelTime[levelTime.size()-level] += omp_get_wtime()-start;
}


void powiter_pprdeg(int &S, int level,int tid, vector<int> &leaf,
                    vector<int> &spprs,vector<int> &spprt,vector<double> &spprv,
                    unsigned int &seed){
//    double start = omp_get_wtime();
    vector<double> pi(graph.n, 0);
    vector<double> residuals(graph.n, 0);
    double rsum=graph.deg[S];
    residuals[S] = graph.deg[S];

    power_iteration(pi, residuals, rsum);

//    fwdTime[fwdTime.size()-level] += omp_get_wtime()-start;
    Fwdidx &fwdidx = fwd_idx[tid];

    int size = leaf.size();

    spprs.reserve(size);
    spprt.reserve(size);
    spprv.reserve(size);

    for(int T:leaf){
        spprs.push_back(S);
        spprt.push_back(T);
        spprv.push_back(pi[T]);
    }
//    nnz=0;
}

void fora_pprdeg(int &S, const string &Father, double rmax, double omega,int level,int tid, vector<bool> &insideLeaf,
                 vector<int> &spprs,vector<int> &spprt,vector<double> &spprv,
                 unsigned int &seed){
//    double start = omp_get_wtime();
    double rsum=0;
    //forward propagation, obtain reserve and residual
    //normal forward push with loading stack
    forward_local_update_linear(S,Father,rsum,rmax, tid);
//    fwdTime[fwdTime.size()-level] += omp_get_wtime()-start;
    Fwdidx &fwdidx = fwd_idx[tid];
    spprs.reserve(fwdidx.first.occur.m_num);
    spprt.reserve(fwdidx.first.occur.m_num);
    spprv.reserve(fwdidx.first.occur.m_num);
//    cout << "fwd time:"<<omp_get_wtime()-start<<endl;
    if (verbose){
//        cout << "fwd time:"<<timeBy(start)<<endl;
        cout <<"rmax: "<<rmax <<endl;
    }

//    start = omp_get_wtime();
    unsigned long long num_random_walk=0;
    compute_ppr_with_fwdidx(S, Father,num_random_walk,omega, rsum, tid,insideLeaf, spprs,spprt,spprv,seed);
//    rwTime[rwTime.size()-level] += omp_get_wtime()-start;
    //    cout<<"walk number: "<<num_random_walk<<endl;
//    cout << "rw time:"<<omp_get_wtime()-start<<endl;
    // clear virtual_fwd_idx
    fwdidx.first.clean();
    fwdidx.second.clean();
    if (verbose){
        cout << "rsum: "<<rsum<<endl;
        cout << "omega: "<<omega<<endl;
        cout << "rsum*omega: "<<rsum*omega<<endl;
//        cout<<"walk time: "<<timeBy(start)<<endl;
//        cout << "nnz: "<<nnz<<endl;
        cout<<"walk number: "<<num_random_walk<<endl;
        cout<<"--------------"<<endl;
    }
//    nnz=0;
}

void taupush_pprdeg(int &S, const string &Father, double rmax, int level, int tid,  vector<bool> &insideLeaf,
                    vector<int> &spprs, vector<int> &spprt, vector<double> &spprv,
                    unsigned int &seed){
//    double start = omp_get_wtime();
    double rsum=0;
    //forward propagation, obtain reserve and residual
    //normal forward push with loading stack
    forward_local_update_linear(S,Father,rsum,rmax, tid);
//    cout<<"ssppr time:"<<omp_get_wtime()-start<<endl;
//    fwdTime[fwdTime.size()-level] += omp_get_wtime()-start;
    Fwdidx &fwdidx = fwd_idx[tid];


    spprs.reserve(fwdidx.first.occur.m_num);
    spprt.reserve(fwdidx.first.occur.m_num);
    spprv.reserve(fwdidx.first.occur.m_num);
//    cout << "fwd time:"<<omp_get_wtime()-start<<endl;
    if (verbose){
//        cout << "fwd time:"<<timeBy(start)<<endl;
        cout <<"rmax: "<<rmax <<endl;
    }
    unordered_map<int,double> pprdegs;
    for(long i=0; i < fwdidx.first.occur.m_num; i++){
        int leaf = fwdidx.first.occur[i];
        double val = fwdidx.first[leaf];
        if (insideLeaf[leaf])
            pprdegs[leaf] += val;
    }
    for(auto pair:pprdegs){
        int T = pair.first;
        spprs.push_back(S);
        spprt.push_back(T);
        spprv.push_back(pair.second);
    }

    fwdidx.first.clean();
    fwdidx.second.clean();
//    nnz=0;
}

void config_pprdeg(int level, vector<int> &leaf, vector<bool> &insideLeaf, double &dmax, double &epsilonR, int &gamma, double &Delta, double &maxval,
                   double &minval, int &nrow) {
    epsilonR= graph.epR;
    gamma= 1;
    Delta= 0;
    maxval= 1 - pow(graph.n, -2 * epsilonR);
    minval= 1 - exp(-2 * epsilonR);
    nrow= leaf.size();
// Delta = \sum_{C\in V_i} d_{sum}(C)/|Leaf(C)|
//    if(visual_mode==INTERACT_MODE){
    int id=0;
    vector<int> &ids = *leaf2id;
    ids.resize(graph.n,-1);
    for (int i:leaf) {
        insideLeaf[i] = true;
        double deg = (double)graph.deg[i];
        Delta += deg;
        dmax = max(dmax,deg);
        ids[i] = id;
        id++;
    }
    // reserve for openmp
//        fwdcache->resize(id);
//    } else{
//        for (int i:graph.nodes) {
//            Delta += (double)graph.deg[i];
//        }
//    }
//    graph.rmaxs=0;
//    graph.rsums=0;
}
void config_pr_pprdeg(int level, vector<int> &leaf, vector<bool> &insideLeaf, double &epsilonR, double &tau, double &maxval,
                      double &minval, int &nrow, vector<int> &large_target) {
    epsilonR= graph.epR;
    double max_tau= tau;
    tau = 0;
    maxval= 1 - pow(graph.n, -2 * epsilonR);
    minval= 1 - exp(-2 * epsilonR);
    nrow= leaf.size();
    int id=0;
    vector<int> &ids = *leaf2id;
    priority_queue<pair<float,int>> sort_tau;
    ids.resize(graph.n,-1);
    for (int i:leaf) {
        insideLeaf[i] = true;
        ids[i] = id;
        if (isFPSN) {
            if (isBPSN) sort_tau.push(pair<float, int>(pr[i], i));
            tau = max(tau,(double)pr[i]);
        }
        else tau = max(tau,sqrt((double)pr[i]));
        id++;
    }
    if (isBPSN){
        while (sort_tau.top().first>max_tau){
            large_target.push_back(sort_tau.top().second);
            sort_tau.pop();
        }
        tau = sort_tau.top().first;
    }
//    cerr<<"threshold tau:"<<max_tau<<endl;
//    cerr<<"local max tau:"<<tau<<endl;
//    cerr<<"prune size:"<<large_target.size()<<endl;

}

void allpair_pprdeg(const string &supernode, int level,
                    vector<int> &leaf, vector<int> &node2id){
    double start,epsilonR;
    double dmax=0;
    int gamma,nrow;
    double Delta,epsilonp,maxval,minval;
    vector<bool> insideLeaf(graph.n, false);
    config_pprdeg(level, leaf, insideLeaf, dmax, epsilonR, gamma, Delta, maxval, minval, nrow);
    start = omp_get_wtime();
    int childsize = leaf.size();
    double total_rinit = 0;
#pragma omp parallel num_threads(thread_num)
#pragma omp for schedule(static) nowait
    for (int id=0;id<leaf.size();id++) {
        double s1 = omp_get_wtime();
        double rmax = 0;
        double omega = 0;
        assert(omp_get_num_threads()==thread_num);
        vector<int> spprs, spprt;
        vector<double> spprv;
        int i = leaf[id];
        // parameter settings for each ss superppr
        // set deltap for current source supernode as the average degree of its leaf nodes
//        double deltap = graph.dbar * graph.delta;
        double deltap = graph.deg[i] * graph.delta;
        // set the corresponding epsilonp for current source supernode
        epsilonp = min(maxval, max(1 - pow(2 * deltap / exp(1), epsilonR), minval));
        if (isFORASN)
            rmax = epsilonp * sqrt( Delta * deltap / (2.0 * graph.m) / (2 + 2 * epsilonp / 3.0) / log(1 / graph.pfail));
        else
            rmax = epsilonp * sqrt(graph.deg[i] * deltap / (2.0 * graph.m) / (2 + 2 * epsilonp / 3.0) / log(1 / graph.pfail));

        omega = (2 + 2 * epsilonp / 3.0) * log(1 / graph.pfail) / epsilonp / epsilonp / deltap ;
//        }
        if (verbose) {
            cout << "Delta: " << Delta << endl;
            cout << "deltap: " << deltap << endl;
            cout << "rinit: "<< graph.g[i].size() << endl;
            cout << "epsilonp: " << epsilonp << endl;
            cout << "deltap*epsilonp: " << deltap * epsilonp << endl;
        }
        unsigned int seed = (id+1)*(omp_get_thread_num()+1);
//        total_rinit += graph.g[i].size();
        fora_pprdeg(i, supernode, rmax, omega, level, omp_get_thread_num(),insideLeaf,spprs,spprt,spprv,seed);

        long leng = spprs.size();
        int idx;
        for (int i = 0; i < leng; ++i) {
            idx = get_indice(node2id[spprs[i]], node2id[spprt[i]], nrow);
            PPRDeg[idx] = spprv[i];
        }
    }
    timeElasped += omp_get_wtime()-start;
//    levelTime[levelTime.size()-level] += omp_get_wtime()-start;
}

void allpair_pprdeg_powiter(const string &supernode, int level,
                            vector<int> &leaf, vector<int> &node2id){
    double start;
    int nrow = leaf.size();
    start = omp_get_wtime();
#pragma omp parallel num_threads(thread_num)
#pragma omp for schedule(static) nowait
    for (int id=0;id<leaf.size();id++) {
        double s1 = omp_get_wtime();
        assert(omp_get_num_threads()==thread_num);
        vector<int> spprs, spprt;
        vector<double> spprv;
        int i = leaf[id];

        unsigned int seed = (id+1)*(omp_get_thread_num()+1);
        powiter_pprdeg(i, level, omp_get_thread_num(),leaf,spprs,spprt,spprv,seed);

        long leng = spprs.size();
        int idx;
        for (int i = 0; i < leng; ++i) {
            idx = get_indice(node2id[spprs[i]], node2id[spprt[i]], nrow);
            PPRDeg[idx] = spprv[i];
        }
    }

    timeElasped += omp_get_wtime()-start;
//    levelTime[levelTime.size()-level] += omp_get_wtime()-start;
}


void allpair_pprdeg_pr(const string &supernode, int level,
                       vector<int> &leaf, vector<int> &node2id){
    double start,epsilonR;
    int nrow;
    double epsilonp,maxval,minval,tau;
    if (isBPSN){
        tau = graph.tau;
    } else{
        tau = 0;
    }

    vector<int> large_target;
    vector<bool> insideLeaf(graph.n, false);
    config_pr_pprdeg(level, leaf, insideLeaf, epsilonR, tau, maxval, minval, nrow, large_target);
    start = omp_get_wtime();
    int childsize = leaf.size();
    double total_rinit = 0;
//    cout<<"leaf: "<<leaf.size()<<endl;
#pragma omp parallel num_threads(thread_num)
#pragma omp for schedule(static) nowait
    for (int id=0;id<leaf.size();id++) {
        assert(omp_get_num_threads()==thread_num);
        vector<int> spprs, spprt;
        vector<double> spprv;
        int i = leaf[id];
        // parameter settings for each ss superppr
        // set deltap for current source supernode as the average degree of its leaf nodes
        double deltap = graph.dbar * graph.delta;
//        double deltap = graph.g[i].size() * graph.delta;
        // set the corresponding epsilonp for current source supernode
        epsilonp = min(maxval, max(1 - pow(2 * deltap / exp(1), epsilonR), minval));

        double rmax;
        rmax= deltap*epsilonp/tau/2/graph.m;
//        cout<<"rmax:"<<rmax<<endl;
        if (verbose) {
            cout << "deltap: " << deltap << endl;
            cout << "rinit: "<< graph.g[i].size() << endl;
            cout << "epsilonp: " << epsilonp << endl;
            cout << "deltap*epsilonp: " << deltap * epsilonp << endl;
        }
        unsigned int seed = (id+1)*(omp_get_thread_num()+1);

        taupush_pprdeg(i, supernode, rmax, level, omp_get_thread_num(), insideLeaf, spprs, spprt, spprv, seed);

        long leng = spprs.size();
        int idx;
        for (int i = 0; i < leng; ++i) {
            idx = get_indice(node2id[spprs[i]], node2id[spprt[i]], nrow);
            PPRDeg[idx] = spprv[i];
        }
        // load the bwdidx here
        for(int t:large_target){
            auto pair = bwd_idx_info[t];
            for (int j = 0; j < pair.second; ++j) {
                int s = leaf[j];
                idx = get_indice(node2id[s], node2id[t], nrow);
                PPRDeg[idx] = bwd_idx_in_cluster[pair.first+j];
            }
        }

    }
    timeElasped += omp_get_wtime()-start;
//    levelTime[levelTime.size()-level] += omp_get_wtime()-start;
}


void store_position(const string &supernode, const Coordinate2d &projections){
    unsigned M = projections.rows();
    string prefix = storepath+supernode+"_"+alg;
    string datapath;
    datapath = prefix+".x";
    ofstream fout(datapath, ios::out | ios::binary);
    fout.write((char*)(projections.data()+0), M * sizeof(double));
    fout.close();
    datapath = prefix+".y";
    ofstream fout1(datapath, ios::out | ios::binary);
    fout1.write((char*)(projections.data()+M), M * sizeof(double));
    fout1.close();
}
void store_radius(const string &supernode, const vector<double> &r){
    unsigned M = r.size();
    string prefix = storepath+supernode+"_"+alg;
    string datapath;
    datapath = prefix+".r";
    ofstream fout(datapath, ios::out | ios::binary);
    fout.write((char*)&r[0], M * sizeof(double));
    fout.close();
}


double calculateSigma( const vector<double>& Weights, const vector<double>& Distances, const Coordinate2d & InitialZ, vector<double>& dists ) {
    unsigned M = InitialZ.rows();
    double sigma=0;
    for(unsigned i=1; i<M; ++i) {
        for(unsigned j=0; j<i; ++j) {
            double dlow=0;
            for(unsigned k=0; k<InitialZ.cols(); ++k) {
                double tmp=InitialZ(i,k) - InitialZ(j,k);
                dlow+=tmp*tmp;
            }
            dists[get_indice(i,j,M)]=dists[get_indice(j,i,M)]=sqrt(dlow);
            double tmp3 = Distances[get_indice(i,j,M)] - dists[get_indice(i,j,M)];
            sigma += Weights[get_indice(i,j,M)]*tmp3*tmp3;
        }
    }
    return sigma;
}


void Smacof( const vector<double>& Weights, const vector<double>& Distances, const double& tol, const unsigned& maxloops, Coordinate2d& InitialZ ) {
    unsigned M = InitialZ.rows();
    if(M>100){
        Eigen::initParallel();
        Eigen::setNbThreads(thread_num);
    }

    // Calculate V
    MatrixXd V(M,M);
    double totalWeight=0.;
    for(unsigned i=0; i<M; ++i) {
        for(unsigned j=0; j<M; ++j) {
            if(i==j) continue;
            V(i,j)=-Weights[get_indice(i,j,M)];
            if( j<i ) totalWeight+=Weights[get_indice(i,j,M)];
        }
        for(unsigned j=0; j<M; ++j) {
            if(i==j)continue;
            V(i,i)-=V(i,j);
        }
    }

    // And pseudo invert V
    MatrixXd mypseudo = pseudoInverse(V,1e-15);
    vector<double> dists(M*M,0);
    double myfirstsig = calculateSigma( Weights, Distances, InitialZ, dists ) / totalWeight;

    // initial sigma is made up of the original distances minus the distances between the projections all squared.
    MatrixXd BZ( M, M );
    Coordinate2d newZ;
    for(unsigned n=0; n<maxloops; ++n) {
//    if(n==maxloops-1) plumed_merror("ran out of steps in SMACOF algorithm");

        // Recompute BZ matrix
        for(unsigned i=0; i<M; ++i) {
            for(unsigned j=0; j<M; ++j) {
                if(i==j) continue;  //skips over the diagonal elements

                if( dists[get_indice(i,j,M)]>0 ) BZ(i,j) = -Weights[get_indice(i,j,M)]*Distances[get_indice(i,j,M)] / dists[get_indice(i,j,M)];
                else BZ(i,j)=0.;
            }
            //the diagonal elements are -off diagonal elements BZ(i,i)-=BZ(i,j)   (Equation 8.25)
            BZ(i,i)=0; //create the space memory for the diagonal elements which are scalars
            for(unsigned j=0; j<M; ++j) {
                if(i==j) continue;
                BZ(i,i)-=BZ(i,j);
            }
        }
        newZ = mypseudo*(BZ*InitialZ);
        //Compute new sigma
        double newsig = calculateSigma( Weights, Distances, newZ, dists ) / totalWeight;
        //Computing whether the algorithm has converged (has the mass of the potato changed
        //when we put it back in the oven!)
        if( fabs( newsig - myfirstsig )<tol ) break;
        myfirstsig=newsig;
        InitialZ = newZ;
    }
}
void add_radius(vector<double> &D, vector<int>& nodeweight, vector<double>& r, double avgsize, double avgdist, double scaler=0.03){
    double scale = scaler*avgdist/sqrt(avgsize);
    double M = nodeweight.size();
    for (int i = 0; i < M; ++i) {
        r[i] = scale*sqrt(nodeweight[i]);
    }
    // update Distance by adding radius
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            D[get_indice(i,j,M)] += r[i]+r[j];
        }
    }

}



void relativeMDS(const double& tol, const unsigned& maxloops, Coordinate2d& projections, vector<double> &r, vector<int>& nodeweight){
    int nrow = projections.rows();
    vector<double> weights(nrow*nrow);
    vector<double> pprdist(nrow*nrow);

    double val;
    double maxval = 2*log2(graph.n);
    double minval = 2;
    double avgsize = 0, avgdist = 0;
    for (int i = 0; i < nrow; ++i) {
        avgsize += nodeweight[i];
        for (int j = 0; j < nrow; ++j) {
            val = PPRDeg[get_indice(i,j,nrow)]+PPRDeg[get_indice(j,i,nrow)];
            if (val==0) val = maxval;
            else val = 1 - log2(val);
            val = min(max(minval, val), maxval);
            avgdist += val;
            pprdist[get_indice(i,j,nrow)] = val;
            weights[get_indice(i,j,nrow)] = 1/val/val;
        }
    }
    avgsize /= nrow;
    avgdist /= (nrow*nrow);
    add_radius(pprdist, nodeweight, r, avgsize, avgdist);
    Smacof( weights, pprdist, tol, maxloops, projections );

}
void zoom_in(const string &supernode, int level, vector<string> &children){

//    num_rw_level=0;
    int nrow;
    vector<int> node2id(graph.n,-1);
    // get the top of fwd cache stack
//    fwdcache = &fwd_stack[fwd_stack_top];
    leaf2id = &leaf2id_stack[fwd_stack_top];
    vector<int> nodeweight;
    if(level==1){
        vector<int> &leaf = super2leaf[supernode];
        nrow = leaf.size();
        PPRDeg.resize(nrow*nrow,0);
        nodeweight.resize(nrow,1);
        for (int i = 0; i < nrow; ++i) {
            node2id[leaf[i]] = i;
        }
        if (isPowerIter)
            allpair_pprdeg_powiter(supernode,level,leaf,node2id);
        else if (isFPSN)
            allpair_pprdeg_pr(supernode,level,leaf,node2id);
        else
            allpair_pprdeg(supernode,level,leaf,node2id);
        if(verbose){
            for(const auto& each:leaf)
                cout<<each<<" ";
            cout<<endl;
        }
    }
    else{
        vector<int> &child = super2super[supernode];
        nrow = child.size();
        PPRDeg.resize(nrow*nrow,0);
        nodeweight.resize(nrow,1);
        for (int i = 0; i < nrow; ++i) {
            node2id[child[i]] = i;
            nodeweight[i] = super2leaf[children[i]].size();
        }
        if (isPowerIter)
            allpair_super_pprdeg_powiter(supernode, level, children, node2id);
        else if (isFPSN)
            allpair_super_pprdeg_pr(supernode,level,children,node2id);
        else
            allpair_super_pprdeg(supernode,level,children,node2id);
//        cout<<"number of children: "<<children.size()<<endl;
        if(verbose){
            for(const auto& each:children)
                cout<<each<<" ";
            cout<<endl;
        }
    }
//    cerr << "ppr time: "<<timeBy(start)<<endl;
    if (embed_on){
        double start2 = omp_get_wtime();
        double tol=1e-9;
        unsigned maxloops=3000;
        Coordinate2d projections = Coordinate2d::Random(nrow,2);
        vector<double> r(nrow);
        relativeMDS(tol,maxloops,projections,r,nodeweight);
        timeElasped += omp_get_wtime()-start2;
        //    cerr << "mds time: "<<omp_get_wtime()-start2<<endl;
        //    cerr << "num of rw: "<<num_rw_level<<endl;
        store_position(supernode,projections);
        // store supernode radius
        store_radius(supernode,r);
    }
    // increase top for stack
    fwd_stack_top++;
    PPRDeg.clear();
}

void full_visualize(){
    vector<int> &leaf = graph.nodes;
    int nrow = leaf.size();
    vector<int> node2id(graph.n);
    std::iota(std::begin(node2id), std::end(node2id), 0);
    vector<int> nodeweight;
    PPRDeg.resize(nrow*nrow,0);
    nodeweight.resize(nrow,1);
    for (int i = 0; i < nrow; ++i) {
        node2id[leaf[i]] = i;
    }

#pragma omp parallel num_threads(thread_num)
#pragma omp for schedule(static) nowait
    for (int id=0;id<leaf.size();id++) {
        double s1 = omp_get_wtime();
        assert(omp_get_num_threads()==thread_num);
        vector<int> spprs, spprt;
        vector<double> spprv;
        int i = leaf[id];
        unsigned int seed = (id+1)*(omp_get_thread_num()+1);
        powiter_pprdeg(i, 1, omp_get_thread_num(),leaf,spprs,spprt,spprv,seed);
        long leng = spprs.size();
        int idx;
        for (int i = 0; i < leng; ++i) {
            idx = get_indice(node2id[spprs[i]], node2id[spprt[i]], nrow);
            PPRDeg[idx] = spprv[i];
        }
    }

    double start2 = omp_get_wtime();
    double tol=1e-3;
    unsigned maxloops=300;
    Coordinate2d projections = Coordinate2d::Random(nrow,2);
    vector<double> r(nrow);
    relativeMDS(tol,maxloops,projections,r,nodeweight);
    cout<<omp_get_wtime()-start2<<endl;
}

void interactive_visualize(vector<string> &path){
    vector<string> childrens;
    int level;
    double start;
//    path = {"c0_l2_11","c0_l1_637"};
//    path = {"c0_l3_0","c0_l2_11","c0_l1_637"};
//    path = {"c0_l5_0","c0_l4_42","c0_l3_1734","c0_l2_11196","c0_l1_889018"};
    path = {"c0_l2_0"};

    for(const string& supernode:path){
//        string supernode;
//        zoom_in("c0_l3_0");
//        cout<<"select the above supernode to zoom in... (type -1 to exit)"<<endl;
//        cin >> supernode;
        parse_children(supernode, level, childrens);
        cerr<<"zoom into supernode "<<supernode<<"..."<<endl;
//        start = omp_get_wtime();
        zoom_in(supernode,level, childrens);
//        levelTime[levelTime.size()-level] += omp_get_wtime()-start;
//        break;
        if (level==1){
//            cerr<<"reach the leaf level: exiting..."<<endl;
            return;
        }
    }
}


#endif //INTERACT_FORA_ALGO_H
