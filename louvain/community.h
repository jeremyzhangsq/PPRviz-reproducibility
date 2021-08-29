#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h> // memcpy
#include <iostream>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <random>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <unistd.h>
#include <array>
#include <queue>
#include <unordered_set>
#include "sparsehash/dense_hash_map"
#include "graph_binary.h"

using namespace std;
using google::dense_hash_map;

const int VectorDefaultSize=20;
template <typename _T>
class iVector
{
public:
    unsigned int m_size;
    _T* m_data;
    unsigned int m_num;

    void free_mem()
    {
        delete[] m_data;
    }

    iVector()
    {
        //printf("%d\n",VectorDefaultSize);
        m_size = VectorDefaultSize;
        m_data = new _T[VectorDefaultSize];
        m_num = 0;
    }
    iVector( unsigned int n )
    {
        if ( n == 0 )
        {
            n = VectorDefaultSize;
        }
//      printf("iVector allocate: %d\n",n);
        m_size = n;
        m_data = new _T[m_size];
        m_num = 0;
    }
    void push_back( _T d )
    {
        if ( m_num == m_size )
        {
            re_allocate( m_size*2 );
        }
        m_data[m_num] = d ;
        m_num++;
    }
    void push_back( const _T* p, unsigned int len )
    {
        while ( m_num + len > m_size )
        {
            re_allocate( m_size*2 );
        }
        memcpy( m_data+m_num, p, sizeof(_T)*len );
        m_num += len;
    }

    void re_allocate( unsigned int size )
    {
        if ( size < m_num )
        {
            return;
        }
        _T* tmp = new _T[size];
        memcpy( tmp, m_data, sizeof(_T)*m_num );
        m_size = size;
        delete[] m_data;
        m_data = tmp;
    }
    void Sort()
    {
        if ( m_num < 20 )
        {
            int k ;
            _T tmp;
            for ( int i = 0 ; i < m_num-1 ; ++i )
            {
                k = i ;
                for ( int j = i+1 ; j < m_num ; ++j )
                    if ( m_data[j] < m_data[k] ) k = j ;
                if ( k != i )
                {
                    tmp = m_data[i];
                    m_data[i] = m_data[k];
                    m_data[k] = tmp;
                }
            }
        }
        else sort( m_data, m_data+m_num );
    }
    void unique()
    {
        if ( m_num == 0 ) return;
        Sort();
        unsigned int j = 0;
        for ( unsigned int i = 0 ; i < m_num ; ++i )
            if ( !(m_data[i] == m_data[j]) )
            {
                ++j;
                if ( j != i ) m_data[j] = m_data[i];
            }
        m_num = j+1;
    }
    int BinarySearch( _T& data )
    {
        for ( int x = 0 , y = m_num-1 ; x <= y ; )
        {
            int p = (x+y)/2;
            if ( m_data[p] == data ) return p;
            if ( m_data[p] < data ) x = p+1;
            else y = p-1;
        }
        return -1;
    }
    void clean()
    {
        m_num = 0;
    }
    void assign( iVector& t )
    {
        m_num = t.m_num;
        m_size = t.m_size;
        delete[] m_data;
        m_data = t.m_data;
    }

    bool remove( _T& x )
    {
        for ( int l = 0 , r = m_num ; l < r ; )
        {
            int m = (l+r)/2;

            if ( m_data[m] == x )
            {
                m_num--;
                if ( m_num > m ) memmove( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                return true;
            }
            else if ( m_data[m] < x ) l = m+1;
            else r = m;
        }
        return false;
    }

    void sorted_insert( _T& x )
    {
        if ( m_num == 0 )
        {
            push_back( x );
            return;
        }

        if ( m_num == m_size ) re_allocate( m_size*2 );

        int l,r;

        for ( l = 0 , r = m_num ; l < r ; )
        {
            int m = (l+r)/2;
            if ( m_data[m] < x ) l = m+1;
            else r = m;
        }

        if ( l < m_num && m_data[l] == x )
        {
            //printf("Insert Duplicate....\n");
            //cout<<x<<endl;
            //      break;
        }
        else
        {
            if ( m_num > l )
            {
                memmove( m_data+l+1, m_data+l, sizeof(_T)*(m_num-l) );
            }
            m_num++;
            m_data[l] = x;
        }
    }

    bool remove_unsorted( _T& x )
    {
        for ( int m = 0 ; m < m_num ; ++m )
        {
            if ( m_data[m] == x )
            {
                m_num--;
                if ( m_num > m ) memcpy( m_data+m, m_data+m+1, sizeof(_T)*(m_num-m) );
                return true;
            }
        }
        return false;
    }

    _T& operator[]( unsigned int i )
    {
        //if ( i < 0 || i >= m_num )
        //{
        //  printf("iVector [] out of range!!!\n");
        //}
        return m_data[i];
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //close range check for [] in iVector if release

};
template <typename _T>
struct iMap
{
    _T* m_data;
    int m_num;
    int cur;
    iVector<int> occur;
    _T nil;
    iMap()
    {
        m_data = NULL;
        m_num = 0;
        //nil = std::make_pair((long)-9,(long)-9);
        //nil = 1073741834;
    }
    iMap(int size){
        initialize(size);
    }
    void free_mem()
    {
        delete[] m_data;
        occur.free_mem();
    }

    void initialize( int n )
    {
        occur.re_allocate(n);
        occur.clean();
        m_num = n;
        nil = -9;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = nil;
        cur = 0;
    }
    void clean()
    {
        for ( int i = 0 ; i < occur.m_num ; ++i )
        {
            m_data[occur[i]] = nil;
        }
        occur.clean();
        cur = 0;
    }

    //init keys 0-n, value as 0
    void init_keys(int n){
        occur.re_allocate(n);
        occur.clean();
        m_num = n;
        nil = -9;
        if ( m_data != NULL )
            delete[] m_data;
        m_data = new _T[m_num];
        for ( int i = 0 ; i < m_num ; ++i ){
            m_data[i] = 0;
            occur.push_back( i );
            cur++;
        }
    }
    //reset all values to be zero
    void reset_zero_values(){
        // for ( int i = 0 ; i < m_num ; ++i )
        // m_data[i] = 0.0;
        memset( m_data, 0.0, m_num*sizeof(_T) );
    }

    void reset_one_values(){
        for ( int i = 0 ; i < m_num ; ++i )
            m_data[i] = 1.0;
        // memset( m_data, 0.0, m_num*sizeof(_T) );
    }

    _T get( int p )
    {
        //if ( p < 0 || p >= m_num )
        //{
        //  printf("iMap get out of range!!!\n");
        //  return -8;
        //}
        return m_data[p];
    }
    _T& operator[](  int p )
    {
        //if ( i < 0 || i >= m_num )
        //{
        //  printf("iVector [] out of range!!!\n");
        //}
        return m_data[p];
    }
    void erase( int p )
    {
        //if ( p < 0 || p >= m_num )
        //{
        //  printf("iMap get out of range!!!\n");
        //}
        m_data[p] = nil;
        cur--;
    }
    bool notexist( int p )
    {
        return m_data[p] == nil ;
    }
    bool exist( int p )
    {
        return !(m_data[p] == nil);
    }
    void insert( int p , _T d )
    {
        //if ( p < 0 || p >= m_num )
        //{
        //  printf("iMap insert out of range!!!\n");
        //}
        if ( m_data[p] == nil )
        {
            occur.push_back( p );
            cur++;
        }
        m_data[p] = d;
    }
    void inc( int p )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("inc some unexisted point\n");
        //}
        m_data[p]++;
    }
    void inc( int p , int x )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("inc some unexisted point\n");
        //}
        m_data[p] += x;
    }
    void dec( int p )
    {
        //if ( m_data[p] == nil )
        //{
        //  printf("dec some unexisted point\n" );
        //}
        m_data[p]--;
    }
    //close range check when release!!!!!!!!!!!!!!!!!!!!
};

void initmaps(int size, int thread_num);
class Community {
public:
    vector<iMap<double>> neigh_map;
    vector<int> neigh_last;
    int thread_num;
    int chunksize;
    Graph g; // network to compute communities for
    int size; // nummber of nodes in the network and size of all vectors

    vector<int> n2c; // community to which each node belongs
    vector<vector<int>> c2n;
    vector<double> in, tot; // used to compute the modularity participation of each community

    // number of pass for one level computation
    // if -1, compute as many pass as needed to increase modularity
    int nb_pass;

    // a new pass is computed if the last one has generated an increase
    // greater than min_modularity
    // if 0. even a minor increase is enough to go for one more pass
    double min_modularity;

    // constructors:
    // reads graph from file using graph constructor
    // type defined the weighted/unweighted status of the graph file
    Community(char *filename, char *filename_w, int type, int nb_pass, double min_modularity, int thread_num);

    // copy graph
    Community(Graph g, int nb_pass, double min_modularity, int thread_num,int chunk);

    // initiliazes the partition with something else than all nodes alone
//    void init_partition(char *filename_part);

    // display the community of each node
    void display();

    // remove the node from its current community with which it has dnodecomm links
    inline void remove(int node, int comm, double dnodecomm);

    // insert the node in comm with which it shares dnodecomm links
    inline void insert(int node, int comm, double dnodecomm);



    // compute the gain of modularity if node where inserted in comm
    // given that node has dnodecomm links to comm.  The formula is:
    // [(In(comm)+2d(node,comm))/2m - ((tot(comm)+deg(node))/2m)^2]-
    // [In(comm)/2m - (tot(comm)/2m)^2 - (deg(node)/2m)^2]
    // where In(comm)    = number of half-links strictly inside comm
    //       Tot(comm)   = number of half-links inside or outside comm (sum(degrees))
    //       d(node,com) = number of links from node to comm
    //       deg(node)   = node degree
    //       m           = number of links
    inline double modularity_gain(int node, int comm, double dnodecomm, double w_degree);
    inline double modularity_gain(int comm, double dnodecomm, double w_degree);
    // compute the set of neighboring communities of node
    // for each community, gives the number of links from node to comm
    void neigh_comm(int node);

    void neigh_comm(int comm_id, vector<int> &community, int tid);
    // compute the modularity of the current partition
    double modularity();

    // displays the graph of communities as computed by one_level
    void partition2graph();

    // displays the current partition (with communities renumbered from 0 to k-1)
    void display_partition();

    // generates the binary graph of communities as computed by one_level
    Graph partition2graph_binary_old(vector<vector<int>> &comm_nodes);
    Graph partition2graph_binary(vector<vector<int>> &comm_nodes);
    Graph partition2graph_binary();
    // compute communities of the graph for one level
    // return true if some nodes have been moved
    bool one_level();

    // compute communities of the graph for one level
    // return true if some nodes have been moved
    void one_level_new(const int &k);
};

inline void
Community::remove(int node, int comm, double dnodecomm) {
    assert(node >= 0 && node < size);
    tot[comm] -= g.weighted_degree(node);
    in[comm] -= 2 * dnodecomm + g.nb_selfloops(node);
    n2c[node] = -1;
}



inline void
Community::insert(int node, int comm, double dnodecomm) {
    assert(node >= 0 && node < size);
    tot[comm] += g.weighted_degree(node);
    in[comm] += 2 * dnodecomm + g.nb_selfloops(node);
    n2c[node] = comm;
}



inline double
Community::modularity_gain(int node, int comm, double dnodecomm, double w_degree) {
    assert(node >= 0 && node < size);

    double totc = (double) tot[comm];
    double degc = (double) w_degree; // the sum of the weightes of the links incident to "node"
    double m2 = (double) g.total_weight;
    double dnc = (double) dnodecomm; // the sum of weights of the links from "node" to nodes in C

    return (dnc - totc * degc / m2);
}

// it is equal to the sum of modularity gain of each single node
// comm: neigbor community
// dnodecomm: the sum of weights of the links from nodes in current comm to nodes in neighbor comm
// w_degree: the sum of the weightes of the links incident to nodes in current comm
inline double
Community::modularity_gain(int comm, double dnodecomm, double w_degree) {

    double totc = (double) tot[comm];
    double degc = (double) w_degree;
    double m2 = (double) g.total_weight;
    double dnc = (double) dnodecomm;

    return (dnc - totc * degc / m2);
}

#endif // COMMUNITY_H
