#pragma once
#include "neighborhood.hpp"
#include "unionFind.hpp"
/*  graph -- a array of edges arranged by
 *  increasing length.  Allows iteration:
 *
 *  init_iteration() places current at
 *    the start of the array
 *  get_next() returns the current edge
 *    and advances to the next one
 */
class graph
{
private:
  //Helper functions
  int existAFalse(bool *visited);
  bool popEdgesGetClusters(list<edge> *ce, list<edge> *poped_edges, int k, cloud *c);
  void pushBackAllEdges(list<edge> *ce, list<edge> *poped_edges);
  std::list<edge> *getAllEdges();
  std::list<std::list<edge>> *getClusters(list<edge> *ce, list<edge> *poped_edges, int best_k, cloud *c);
  std::list<std::list<edge>> *MSDR(cloud *c, int max_k, std::list<edge> *ce);
  std::list<std::list<edge>> *getBestk(cloud *c, int max_k, std::list<edge> *ce);

public:
  cloud *c;
  bool **isEdge;
  int size, numEdges;
  neighborhood *Adjlist;
  graph(int _size, ifstream *inFile, ifstream *weightFile);
  graph(cloud &_c);
  ~graph()
  {
    delete[] Adjlist;
    delete[] isEdge;
  }

  //Task 1 and 2
  std::list<edge> Boruvska();
  std::list<edge> Prim();
  std::list<edge> *Kruskal();
  //End Task 1 and 2

  //Task 5
  std::list<std::list<edge>> *kClustering(int k, cloud *c);
  //End Task 5
  
  //Task 6
  std::list<std::list<edge>> *getClusters(list<edge> *ce);
  std::list<std::list<edge>> *getBestk(cloud *c, int max_k);
  std::list<std::list<edge>> *MSDR(cloud *c, int max_k);
  //End Task 6

  static double getTotalWeight(std::list<edge> l);

  void refineClustering(cloud *c, list<edge> *ce);
  void graphToFile(std::ofstream *outfile);
  void setLabel(std::list<std::list<edge>> *ret, cloud *c);
  int computeMinima(double *y, int N);
  int get_size()
  {
    return size;
  }
};
