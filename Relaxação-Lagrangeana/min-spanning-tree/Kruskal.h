#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>

using namespace std;

typedef pair<int, int> ii;
typedef vector <vector<double> >vvi;
typedef vector<ii> vii;

class Kruskal{
public:
	Kruskal(vvi dist);

	double MST(size_t nodes);
	vii getEdges();
	vii edges;


private:
	priority_queue <pair<double,ii> > graph;
	vector <size_t> pset;

	void initDisjoint(size_t n);
	size_t findSet(size_t i);
	void unionSet(size_t i, size_t j);
	bool isSameSet(size_t i, size_t j);
};

#endif