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
#define KMAX 10
#define MAXCOST 99999999
#define INITEPS 0.5

typedef pair<int, int> ii;
typedef vector<vector<double>> vvi;
typedef vector<ii> vii;

class Node{
public:
	Node(){};
	Node(double lb, double ub, const vector<double>& _lambda, vector<double>& _curr_lambda, vector<pair<size_t,size_t>>& _best_edges, vector<pair<size_t,size_t>>& _curr_edges);
	~Node();

	void Solve(const vector<vector<double>>& og_cost, vector<vector<double>>& cost, vector <size_t>& _pset, vector<ii>& _non_zero_edges);
	
	double UB;
	vector<pair<int, int>> forbidden_edges = vector<pair<int, int>>();
	
	size_t biggest = 0;
	bool feasible = false;
	bool improved = false;
	double LB;
	vector<double> lambda;
	vector<pair<size_t,size_t>>* best_edges;
	// string code = "";
	
private:
	double real_cost = 0;
	vector<double>* curr_lambda;
	// double curr_LB;
	vector<pair<size_t,size_t>>* curr_edges;

	// vector<ii>* non_zero_edges;
	// vector <size_t>* pset;
	// const vvi* dist;
};

#endif