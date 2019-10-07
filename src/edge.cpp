#include <iostream>
#include <list>
#include <iterator>
#include <math.h>
#include "edge.hpp"
/* edge -- method implementations */
cloud *edge::c = NULL;
bool edge::compare(edge e1, edge e2)
{
  if (e1.weight < e2.weight)
    return true;
  else if (e1.weight == e2.weight)
  {
    if (e1.s < e2.s)
      return true;
    else if (e1.s == e2.s)
    {
      if (e1.t < e2.t)
        return true;
    }
  }
  return false;
}

void edge::print(std::list<edge> l)
{
  std::list<edge>::iterator it;
  for (it = l.begin(); it != l.end(); it++)
  {
    std::cout << "(" << it->get_source() << ", " << it->get_target() << ") ";
  }
  std::cout << '\n';
}

bool edge::set_cloud(cloud *_c)
{
  if (c != NULL)
    return false;

  c = _c;
  return true;
}
double edge::mean(std::list<edge> T){
  std::list<edge>::iterator it;
  double m = 0;
  int count = 0;
  for(it = T.begin(); it != T.end(); it++){
    m+=it->get_weight();
    count++;
  }
  return m/count;
}
double edge::std_deviation(std::list<edge> T){
  std::list<edge>::iterator it;
  double s = 0;
  double m = mean(T);
  int count = 0;
  for(it = T.begin(); it != T.end(); it++){
    s+=(it->get_weight()-m)*(it->get_weight()-m);
    count++;
  }
  return sqrt(s/count);
}
double edge::std_deviation(std::list<std::list<edge>> F){
  //It's what we want to minimize
  std::list<std::list<edge>>::iterator it;
  double s_F = 0;
  double TW = 0;
  for(it = F.begin(); it != F.end(); it++){
    s_F += mean(*it)*(it->size()*std_deviation(*it));
    TW +=mean(*it)*(it->size());
  }
  return s_F/TW;
}

bool edge::contains(std::list<edge> *removed_edges, edge beRemoved){
  std::list<edge>::iterator it;
  for(it = removed_edges->begin(); it != removed_edges->end(); it++){
    if(it->s == beRemoved.s && it->t == beRemoved.t){
      return true;
    }
  }
  return false;
}
