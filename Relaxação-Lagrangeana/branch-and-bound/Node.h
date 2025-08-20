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
#define MINEPS 1e-6
#define KMAX 15
#define MAXCOST 99999999
#define INITEPS 0.5

typedef pair<int, int> ii;
typedef vector<vector<double>> vvi;
typedef vector<ii> vii;

class Node{
public:
	Node(){};
	Node(double lb, double ub, const vector<double>& _lambda/*, const vector<pair<size_t,size_t>>& _edges*/, const string& _code) { LB = lb; UB = ub; curr_lambda = _lambda; /*edges = _edges;*/ code = _code; };

	void Solve(const vector<vector<double>>& og_cost, vector<vector<double>>& cost);
	
	double UB;
	vector<pair<int, int>> forbidden_edges = vector<pair<int, int>>();
	
	size_t biggest = 0;
	bool feasible = false;
	double LB;
	vector<double> lambda;
	vector<pair<size_t,size_t>> edges = vector<pair<size_t,size_t>>();
	string code = "";
	
private:
	double real_cost = 0;
	vector<double> curr_lambda;
	// double curr_LB;
	vector<pair<size_t,size_t>> curr_edges;
};

#endif