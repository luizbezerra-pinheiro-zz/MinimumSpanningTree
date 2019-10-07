#pragma once
#include <list> 
#include <iterator>
#include "cloud.hpp"
#include "cassert"
/* edge -- pairs of cloud point indices with methods to
 *  - compute the distance from source to target
 *  - compare the lengths of two edges -- needed for sorting
 */

class edge {
  int  s, t;

  static cloud *c;  
public:
  double weight; 
  bool operator == (const edge& that) const { return (s == that.s && t == that.t) || (s == that.t && t == that.s); }
  bool operator != (const edge& that) const { return !operator==(that); }
    


  edge(int _s, int _t, int _weight) {
    //assert(c != NULL);
    //assert( _s >= 0 && _s < c->get_n() );
    //assert( _t >= 0 && _t < c->get_n() );
    assert(_weight >=0);
    s = _s;
    t = _t;
    weight = _weight;
  }
  edge(int _s, int _t){
    assert(c != NULL);
    assert( _s >= 0 && _s < c->get_n() );
    assert( _t >= 0 && _t < c->get_n() );
    s = _s;
    t = _t;
    weight = c->get_point(s).dist(c->get_point(t));
  }
  edge(){
    int s=0, t=0;
  }
  //getters
  int get_source() { return s; }
  int get_target() { return t; }
  int get_source() const { return s; }
  int get_target() const { return t; }
  double get_weight(){return weight;}
  double get_weight() const {return weight;}

  //setters
  
  static bool set_cloud (cloud *_c);

  //helper functions
  static double mean(std::list<edge> T);
  static double std_deviation(std::list<edge> T);
  static double std_deviation(std::list<std::list<edge>> F);
  static bool compare (edge e1, edge e2);
  static void print(std::list<edge> l);
  static bool contains(std::list<edge> * removed_edges, edge beRemoved);
};
struct Comp{
    bool operator()(const edge& a, const edge& b){
        return a.get_weight()>b.get_weight();
    }
};
