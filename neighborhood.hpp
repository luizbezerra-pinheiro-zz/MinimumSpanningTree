#pragma once
#include <list> 
#include <iterator>
#include "edge.hpp"

using namespace std; 
/*  DESCRIÇÃO  */
class neighborhood {
  public: 
  list <edge> neighbors;
    neighborhood();
  neighborhood(list <edge> *l);
  ~neighborhood();
  void addEdge(edge * _e);
  void print(); //for testing
};