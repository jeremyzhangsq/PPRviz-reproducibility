#include "lib.h"

void parameter(int argc, char** argv, unordered_map<string, string> &map){
    for(int i=1;i<argc;i+=2){
        map.insert({argv[i], argv[i+1]});
    }
}

void param_config(string &alg){
    if(alg=="powiter"){ isPowerIter =1;}
    else if(alg=="fora"){ isRWIdx = 1;}
    else if(alg=="foraplus"){ isRWIdx=1; isOptFORA = 1;}
    else if(alg=="foraresacc"){ isResAcc = 1;isOnline= true;}
    else if(alg=="foratp"){ isRWIdx=1; isFORASN = 1;}
    else if(alg == "fpsn"){ isFPSN=1;isBPSN=0;}
    else if(alg == "taupush"){ isFPSN=1;isBPSN=1;}
    else{exit(-1);}
}
int main(int argc,char *argv[]) {
    unordered_map<string, string> param;
    parameter(argc, argv, param);
    string filename;
    int buildflag;
    int sample;
    double alpha;
    int random_query;
    int k;
    filename = param.count("-f")?param["-f"]:"amazon";
    buildflag = stoi(param.count("-build")?param["-build"]:"0");
    verbose = stoi(param.count("-verbose")?param["-verbose"]:"0");
    alpha = stof(param.count("-a")?param["-a"]:"0.2");
    k = stoi(param.count("-k")?param["-k"]:"25");
    sample = stoi(param.count("-sample")?param["-sample"]:"1");
    thread_num = 1; // multi-threading is not support currently
    alg = param.count("-alg")?param["-alg"]:"taupush";
    random_query = stoi(param.count("-random")?param["-random"]:"1");
    embed_on = stoi(param.count("-embed")?param["-embed"]:"0");
    param_config(alg);
    int seed = stoi(param.count("-seed")?param["-seed"]:"2");
    srand(seed);
    string datapath = "../dataset/"+filename;

//    if (full_mode){
//        vector<int> ks = {25,50,100,500,1000};
//        for (int k:ks){
//            datapath = "../dataset/rand_"+to_string(k);
//            graph = Graph(datapath,alpha);
//            full_visualize();
//        }
//        return 0;
//    }

    if (verbose)
        cout << "dataset: "<<filename<<endl;
    hiename = "../louvain/hierachy-output/"+filename +"_"+to_string(k)+".dat";
    mapname = "../louvain/mapping-output/"+filename +"_"+to_string(k)+".dat";
    rootname = "../louvain/hierachy-output/"+filename +"_"+to_string(k)+".root";
    storepath = "../"+filename+"_idx/"+filename+"ds250" +"_"+to_string(k);

    graph = Graph(datapath,alpha,k);
    int max_level = load_multilevel();
    graph.max_level = max_level;

//#ifdef linux
//    cerr << "Graph memory usage (MB):" << getMemory()<< endl << endl;
//#endif

    prpath = "../pr_idx/"+filename+".dnpr";
    if ((!buildflag && isFPSN) or (isBPSN)){
        prpath = "../pr_idx/"+filename+".dnpr";
        deserialize_pr();
    }

    if (buildflag){
        if (isBPSN)
            rwpath = "../bwd_idx/"+filename;
        if (isRWIdx)
            rwpath = "../rwidx/"+filename+"randwalks";

//        int threads[] = {64,32,16,8,4,2,1};
        int threads[] = {1};
        for(int each:threads) {
            thread_num = each < omp_get_max_threads() ? each : omp_get_max_threads();

            if (isBPSN)
                build_bwdpush();
            if (isFPSN and !isBPSN)
                build_dnpr();
            if (isRWIdx)
                build_rwidx();
        }
    }
    else{
        // load rwidx
        if (isRWIdx){
            rwpath = "../rwidx/"+filename+"randwalks"+"_"+alg;
            deserialize_idx();
        }
        if (isBPSN){
            rwpath = "../bwd_idx/"+filename;
            deserialize_bwd();
        }
        if (!random_query){
            if (!isFPSN){
                prpath = "../pr_idx/"+filename+".dnpr";
                deserialize_pr();
            }
            top_k_hub_cluster(sample);
        }

        int threads[] = {1};
        for(int each:threads){
            thread_num = each<omp_get_max_threads()? each:omp_get_max_threads();
            double totaldnpprtime = 0;
            double totalembedtime = 0;
            for (int i = 0; i < sample; ++i) {
                timeElasped = 0;
                embedTimeElapsed = 0;
                vector<string> path;
                init_container();
                if (random_query)
                    generate_random_path(path,max_level);
                else{
                    path.emplace_back(hubcluster[i]);
                }
                interactive_visualize(path);
                cerr <<timeElasped<<endl;
//                cout<<(timeElasped-embedTimeElapsed)<<endl;
                totaldnpprtime += (timeElasped-embedTimeElapsed);
                totalembedtime += embedTimeElapsed;
                if (!random_query){
                    int size = super2leaf[hubcluster[i]].size();
                    cout<<timeElasped/size<<endl;
                }
            }
            if (random_query)
                cout<<totaldnpprtime/sample<<" "<<totalembedtime/sample<<endl;
        }

    }
//#ifdef linux
//    cerr << "Total memory usage (MB):" << getMemory()<< endl << endl;
//#endif
    return 0;
}
