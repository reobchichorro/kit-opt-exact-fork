#ifndef NODE_H
#define NODE_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <utility>
#include <numeric>
#include <cmath>
#include "Kruskal.h"

using namespace std;

#define EPS 1e-6
#define MINEPS 1e-5
#define KMAX 30
#define MAXCOST 99999999
#define INITEPS 1

typedef pair<int, int> ii;
typedef vector<vector<long double>> vvi;
typedef vector<ii> vii;

class Node{
public:
	Node(){};
	Node(long double lb, long double ub, const vector<long double>& _lambda, vector<long double>& _curr_lambda, vector<pair<size_t,size_t>>& _best_edges, vector<pair<size_t,size_t>>& _curr_edges);
	~Node();

	void Solve(const vector<vector<long double>>& og_cost, vector<vector<long double>>& cost, vector <size_t>& _pset, vector<ii>& _non_zero_edges, vector<int>& g);
	
	long double UB;
	vector<pair<int, int>> forbidden_edges = vector<pair<int, int>>();
	
	size_t biggest = 0;
	bool feasible = false;
	bool improved = false;
	long double LB;
	vector<long double> lambda;
	vector<pair<size_t,size_t>>* best_edges;
	// string code = "";
	
private:
	long double real_cost = 0;
	vector<long double>* curr_lambda;
	// long double curr_LB;
	vector<pair<size_t,size_t>>* curr_edges;

	// vector<ii>* non_zero_edges;
	// vector <size_t>* pset;
	// const vvi* dist;
};

#endif