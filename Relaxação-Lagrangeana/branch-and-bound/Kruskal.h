#ifndef KRUSKAL_H
#define KRUSKAL_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <queue>
#include <utility>
#include <algorithm>

using namespace std;

typedef pair<int, int> ii;
typedef vector<vector<double>> vvi;
typedef vector<ii> vii;

class Kruskal{
public:
	Kruskal(const vvi& _dist, int zero_setup, vector <size_t>& _pset, vector<ii>& _non_zero_edges, vector<pair<size_t,size_t>>* _edges);
	~Kruskal();

	double MST(size_t nodes);
	vector<pair<size_t,size_t>>* edges;

private:
	vector<ii>* non_zero_edges;
	vector <size_t>* pset;
	const vvi* dist;

	void initDisjoint(size_t n);
	size_t findSet(size_t i);
	void unionSet(size_t i, size_t j);
	bool isSameSet(size_t i, size_t j);
};

#endif