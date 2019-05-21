#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "cloud.hpp"

cloud::cloud(int _d, int _nmax, int _k)
{
  d = _d;
  point::d = _d;
  n = 0;
  k = _k;

  nmax = _nmax;

  points = new point[nmax];
  centers = NULL;
}
cloud::~cloud()
{
  delete[] points;
  delete[] centers;
}

void cloud::add_point(point &p, int label)
{
  assert(n < nmax);

  for (int m = 0; m < point::get_dim(); m++)
  {
    points[n].coords[m] = p.coords[m];
  }
  points[n].label = label;
  n++;
}

point &cloud::get_point(int i)
{
  return points[i];
}

void cloud::set_center(point &p, int j)
{
  for (int m = 0; m < d; m++)
    centers[j].coords[m] = p.coords[m];
}

double cloud::intracluster_variance()
{
  double sum = 0.0;
  for (int i = 0; i < n; i++)
  {
    if (points[i].label != 0)
      sum += points[i].dist(centers[points[i].label]) *
             points[i].dist(centers[points[i].label]);
  }

  return sum / n;
}
double cloud::intracluster_distance_m()
{
  double sum = 0.0;
  for (int i = 0; i < n; i++)
  {
    if (points[i].label != 0)
      sum += points[i].dist(centers[points[i].label]);
  }
  return sum / n;
}
double cloud::intracluster_distance()
{
  double max = __DBL_MIN__;
  for (int i = 0; i < n; i++)
  {
    for (int j = i + 1; j < n; j++)
    {
      if (points[i].label != 0 && points[j].label != 0 && (points[i].label == points[j].label))
      {
        if (points[i].dist(points[j]) > max)
        {
          max = points[i].dist(points[j]);
        }
      }
    }
  }
  if (max == __DBL_MIN__)
    return (__DBL_MAX__);
  return max;
}
double cloud::intercluster_distance()
{
  double min = __DBL_MAX__;
  for (int i = 1; i <= k; i++) // <= because centers go from [0..k] inclus because centers[0] keeps the trash
  {
    for (int j = i + 1; j <= k; j++)
    {
      if (centers[i].dist(centers[j]) < min)
        min = centers[i].dist(centers[j]);
    }
  }
  if (min != __DBL_MAX__)
    return min;
  else if (k == 1)
  {
    return 0.1; //It's not that bad case
  }
  else
    return (__DBL_MIN__);
}
double cloud::intercluster_distance_m()
{
  double m = 0;
  for (int i = 1; i < k; i++) // <= because centers go from [0..k] inclus because centers[0] keeps the trash
  {
    double min = __DBL_MAX__;
    double ant = min;
    for (int j = i + 1; j <= k; j++)
    {
      if(centers[i].dist(centers[j]) < min){
        ant = min;
        min = centers[i].dist(centers[j]);
      }
    }
    if(ant == __DBL_MAX__) ant = min;
    m+=ant;
  }
  return m/k;
  
}

double cloud::intercluster_std(){
  double m = intercluster_distance_m();
  double s =0;
  for (int i = 1; i < k; i++) // <= because centers go from [0..k] inclus because centers[0] keeps the trash
  {
    double min = __DBL_MAX__;
    double ant = min;
    for (int j = i + 1; j <= k; j++)
    {
      if(centers[i].dist(centers[j]) < min){
        ant = min;
        min = centers[i].dist(centers[j]);
      }
    }
    if(ant = __DBL_MAX__) ant = min;
    s+=(ant-m)*(ant-m);
  }
  s = s/k;
  return sqrt(s);
}

void cloud::set_centroid_centers()
{
  int *counts = new int[k + 1]; //counts[0] just count how many ignored points we have
  for (int j = 0; j <= k; j++)
    counts[j] = 0;
  for (int i = 0; i < n; i++)
  {
    counts[points[i].label]++;
  }

  for (int j = 0; j <= k; j++)
    if (counts[j] != 0)
      for (int m = 0; m < d; m++)
        centers[j].coords[m] = 0.0;

  for (int i = 0; i < n; i++)
  {
    for (int m = 0; m < d; m++)
    {
      centers[points[i].label].coords[m] += points[i].coords[m];
    }
  }

  for (int j = 0; j <= k; j++)
    if (counts[j] != 0)
      for (int m = 0; m < d; m++)
        centers[j].coords[m] /= counts[j];

  delete[] counts;
}

double cloud::utilityFunction()
{
  set_centroid_centers();
  double intrac = intracluster_distance();
  double interc = intercluster_distance();
  if (intrac == __DBL_MAX__ || interc == __DBL_MIN__)
  {
    return __DBL_MAX__;
  }
  return intracluster_distance() / (intercluster_distance()*k*k);
}

void cloud::load(std::ifstream &is)
{
  std::string line;

  assert(is.is_open());

  // point to read into
  point p;
  std::string trash;
  // labels to cycle through
  int label = 0;
  while (getline(is, line))
  {
    std::stringstream s(line);
    std::string word;
    int _s, _t;
    // read new points
    for (int m = 0; m < d; m++)
    {
      s >> p.coords[m];
    }

    add_point(p, label);
  }
  k = 0;
  is.close();
}

void cloud::save(std::ofstream &out)
{
  assert(out.is_open());
  
  if(k==0)
  {
    out << "Data not analysed. \n";
  }
  for (int i = 0; i < n; i++)
  {
    for (int m = 0; m < d; m++)
    {
      out << points[i].coords[m] << " ";
    }
    out << points[i].label << '\n';
  }
  out.close();
}