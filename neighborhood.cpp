#include <iostream>
#include "neighborhood.hpp"

neighborhood::neighborhood()
{
    neighbors = list<edge>();
};
neighborhood::neighborhood(list <edge> *l){
    neighbors = *l;
};
neighborhood::~neighborhood()
{
    neighbors.clear();
};

void neighborhood::addEdge(edge *_e)
{
    neighbors.push_back(*_e);
};
void neighborhood::print()
{
    //Create an iterator of std::list
    std::list<edge>::iterator it;
    // Make iterate point to begining and incerement it one by one till it reaches the end of list.
    it = neighbors.begin();
    std::cout << it->get_source() << " edges = [" << it->get_target();
    it++;
    for (; it != neighbors.end(); it++)
    {
        // Access the object through iterator
        int t = it->get_target();
        //Print the contents
        std::cout << ", " << t;
    }
    std::cout << "]\n";
};