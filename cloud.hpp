#pragma once // single incl.

#include <fstream>
#include <iostream>
#include <list>
#include "point.hpp"

class cloud
{
  int n;
  int d;
  int k;

  // maximum possible number of points
  int nmax;

  point *points;
  point *centers;
  int *representants;

public:
  cloud(int _d, int _nmax, int _k);
  ~cloud();

  // Getters

  int get_n() { return n; }
  int get_k() { return k; }
  int get_d() { return d; }
  int get_representant(int label)
  {
    if (label <= k)
      return representants[label];
    else
      return 0;
  };
  point &get_point(int i);
  point &get_center(int j)
  {
    if (j <= k)
      return centers[j];
  };

  //Setters
  void set_center(point &p, int j);
  void set_k(int _k)
  {
    k = _k;
    delete[] centers;
    centers = new point[k + 1]; //the label 0 is to keep the trash.
  }
  void set_label()
  {
    for (int i = 0; i < n; i++)
    {
      points[i].label = 0;
    }
  }
  void set_representants()
  {
    representants = new int[k + 1];
    for(int i=0; i<=k; i++){
      representants[i] = 0;
    }
    for (int i = 0; i < n; i++)
    {
      for (int r = 1; r <= k; r++)
      {
        if (points[representants[r]].dist(points[i]) < points[representants[r]].dist(centers[r]))
        {
          representants[r] = i;
        }
      }
    }
  }
  void set_centroid_centers();

  // Helper methods

  void add_point(point &p, int label);
  double intracluster_variance();
  double intracluster_distance_m();
  double intracluster_distance();

  double intercluster_std();
  double intercluster_distance();
  double intercluster_distance_m();
  double utilityFunction();

  void load(std::ifstream &is);
  void save(std::ofstream &out);

  void print_labels()
  {
    for (int i = 0; i < n; i++)
    {
      if (points[i].label != 0)
        std::cout << points[i].label << "\n";
    }
  }
};