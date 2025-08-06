#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>

using namespace std;

class Node {
public:
	Node(){};
	double lower_bound; // cost of MST
	double real_cost = 0;
	vector<pair<int, int>> forbidden_edges;
	vector<double> lambda;
	size_t biggest = 0; // vtx with biggest degree
	double epsilon = 0.125;
	size_t k = 0;

private:
};

#endif