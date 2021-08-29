#ifndef GRAPH_H
#define GRAPH_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include "omp.h"
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

class Graph {
 public:
  int nb_nodes;
  long nb_links;
  double total_weight;  

  vector<long> degrees;
  vector<int> links;
  vector<float> weights;
    vector<int> old2new;
    vector<int> new2old;
  Graph();

  Graph(char *filename, char *filename_w, int type);
  
  Graph(int nb_nodes, int nb_links, double total_weight, int *degrees, int *links, float *weights);

  Graph(Graph& orginG, vector<int> &nodes);

  void display(void);
  void display_reverse(void);
  void display_binary(char *outfile);
  bool check_symmetry();


  // return the number of neighbors (degree) of the node
  inline int nb_neighbors(int node);

  // return the number of self loops of the node
  inline double nb_selfloops(int node);

  // return the weighted degree of the node
  inline double weighted_degree(int node);

    // return the weighted degree of current community
  inline double weighted_degree(vector<int> &community);

  // return pointers to the first neighbor and first weight of the node
  inline pair<vector<int>::iterator, vector<float>::iterator > neighbors(int node);
};


inline int
Graph::nb_neighbors(int node) {
  assert(node>=0 && node<nb_nodes);

  if (node==0)
    return degrees[0];
  else
    return degrees[node]-degrees[node-1];
}

inline double
Graph::nb_selfloops(int node) {
  assert(node>=0 && node<nb_nodes);

  pair<vector<int>::iterator, vector<float>::iterator > p = neighbors(node);
  for (int i=0 ; i<nb_neighbors(node) ; i++) {
    if (*(p.first+i)==node) {
      if (weights.size()!=0)
	return (double)*(p.second+i);
      else 
	return 1.;
    }
  }
  return 0.;
}

inline double
Graph::weighted_degree(int node) {
  assert(node>=0 && node<nb_nodes);

  if (weights.size()==0)
    return (double)nb_neighbors(node);
  else {
    pair<vector<int>::iterator, vector<float>::iterator > p = neighbors(node);
    double res = 0;
    for (int i=0 ; i<nb_neighbors(node) ; i++) {
      res += (double)*(p.second+i);
    }
    return res;
  }
}

inline double
Graph::weighted_degree(vector<int> &community) {
    double weight = 0;
    for(int node:community)
        weight += weighted_degree(node);
    return weight;
}

inline pair<vector<int>::iterator, vector<float>::iterator >
Graph::neighbors(int node) {
  assert(node>=0 && node<nb_nodes);

  if (node==0)
    return make_pair(links.begin(), weights.begin());
  else if (weights.size()!=0)
    return make_pair(links.begin()+degrees[node-1], weights.begin()+degrees[node-1]);
  else
    return make_pair(links.begin()+degrees[node-1], weights.begin());
}


#endif // GRAPH_H
