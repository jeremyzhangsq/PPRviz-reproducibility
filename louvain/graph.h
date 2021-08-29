#ifndef GRAPH_H
#define GRAPH_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#define WEIGHTED   0
#define UNWEIGHTED 1

using namespace std;

class Graph {
 public:
  vector<vector<pair<int,float> > > links;
  
  Graph (char *filename, int type);
  
  void clean(int type);
  void renumber(int type);
  void display(int type);
  void display_binary(char *filename, char *filename_w, int type);
};

#endif // GRAPH_H
