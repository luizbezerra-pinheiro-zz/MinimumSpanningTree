#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <chrono>
#include "graph.hpp"
#include "neighborhood.hpp"
#include "edge.hpp"
using namespace std::chrono;

int main(int argc, char *argv[])
{
    const char *infile = argv[1];
    const char *outfile = argv[2];
    ifstream input(infile);
    ifstream weightFile;
    ofstream out;
    weightFile.open("input/weightFile1.txt");
    if (!input)
    {
        cerr << "Unable to open file data file \n";
        exit(1); // call system to stop
    }
    if (!weightFile)
    {
        cerr << "Unable to open file weight file \n";
        exit(1); // call system to stop
    }
    int N;
    string line;
    for (int i = 0; i < 17; i++)
    {
        //read the header
        getline(input, line);
        if (i == 4)
        {
            std::stringstream s(line);
            std::string word;
            s >> word;
            s >> word;
            s >> word;
            std::stringstream(word) >> N;
        }
    }
    std::cout << N << " oi\n";
    auto t1 = high_resolution_clock::now();
    graph *a = new graph(N, &input, &weightFile);
    auto t2 = high_resolution_clock::now();
    auto durationInitGraph = duration_cast<microseconds>(t2 - t1);
    std::list<edge> kruskal = *(a->Kruskal());
     t1 = high_resolution_clock::now();
    auto durationKruskal = duration_cast<microseconds>(t1 - t2);
    std::list<edge> prim = a->Prim();
     t2 = high_resolution_clock::now();
    auto durationPrim = duration_cast<microseconds>(t2 - t1);
    std::list<edge> borusvska = a->Boruvska();
     t1 = high_resolution_clock::now();
    auto durationBoruvska = duration_cast<microseconds>(t1 - t2);
    std::cout << "Weight MST BORUVSKA : " << graph::getTotalWeight(borusvska) << " In: " << durationBoruvska.count() << "us\n";
    std::cout << "Weight MST KRUSKAL : " << graph::getTotalWeight(kruskal) << " In: " << durationKruskal.count()<< "us\n";
    std::cout << "Weight MST PRIM : " << graph::getTotalWeight(prim)<< " In: " << durationPrim.count() << "us\n";

    return 0;
}