#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>
#include "Node.h"

using namespace std;

typedef pair<int, int> ii;
typedef vector<vector<double>> vvi;
typedef vector<ii> vii;

#define USE_ZERO 0
#define TWO_CLOSEST 1

class Kruskal{
public:
	Kruskal(const vvi& dist, int zero_setup);

	double MST(size_t nodes);
	vector<pair<size_t,size_t>> getEdges();
	vector<pair<size_t,size_t>> edges;

	int zero_setup;

private:
	// priority_queue <pair<double,ii> > graph;
	vector<pair<double, ii>> non_zero_edges;
	vector <size_t> pset;
	vector<double> zero_edges;

	void initDisjoint(size_t n);
	size_t findSet(size_t i);
	void unionSet(size_t i, size_t j);
	bool isSameSet(size_t i, size_t j);
};

#endif