#pragma once

class point
{
public:
	static int d;
	double *coords;
	int label;

	point();

	~point();

	static int get_dim ();
	
	static bool set_dim (int _d);

	void print();

	double dist(point &q);
};

