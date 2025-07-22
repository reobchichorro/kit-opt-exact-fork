#include <iostream>
using namespace std;

#include "data.h"
#include "hungarian.h"
#include <vector>
#include <list>
#include <unordered_map>

struct Node {
	vector<pair<int, int>> forbidden_arcs = vector<pair<int,int>>();
	vector<vector<int>> subtours = vector<vector<int>>();
	double lower_bound; // cost of hungarian solution
	int chosen = -1; // index of the smallest subtours
	bool feasible = true;
};

void findSubtours(hungarian_problem_t *p, Data *data, Node* node) {
	vector<int> disjointedSets = vector<int>(data->getDimension());

	for (size_t i = 0; i < disjointedSets.size(); i++)
		disjointedSets[i] = i;
	for (size_t i = 0; i < disjointedSets.size(); i++) {
		disjointedSets[p->successors[i]] = disjointedSets[i]; // TODO: double-check correctness of this
	}
	// for (size_t i = 0; i < disjointedSets.size(); i++) { // TODO: this can be done in O(n) instead of O(n^2)
	// 	for (size_t j = 0; j < disjointedSets.size(); j++) {
	// 		if (p->assignment[i][j] == 1) {
	// 			disjointedSets[j] = disjointedSets[i]; // TODO: double-check correctness of this
	// 			successors[i] = j;
	// 		}
	// 	}
	// }
	unordered_map<int, size_t> subtoursSizes;
	for (size_t i = 0; i < disjointedSets.size(); i++) {
		if (subtoursSizes.contains(disjointedSets[i]))
			subtoursSizes[disjointedSets[i]]++;
		else
			subtoursSizes.emplace(disjointedSets[i], 1);
	}
	size_t smallest = disjointedSets.size() + 1;
	for (auto it = subtoursSizes.cbegin(); it != subtoursSizes.cend(); it++) {
		if ((*it).second < smallest) {
			node->chosen = node->subtours.size();
			smallest = (*it).second;
		}
		node->subtours.push_back(vector<int>());
		int idx = (*it).first;
		do {
			node->subtours.back().push_back(idx);
			idx = p->successors[idx];
		} while(idx != (*it).first);
	}
	if (smallest < disjointedSets.size())
		node->feasible = false;
	hungarian_free(p);
}

void hungarian(hungarian_problem_t *p, Data *data, double **cost, Node* node) {
	hungarian_init(p, cost, data->getDimension(), data->getDimension(), HUNGARIAN_MODE_MINIMIZE_COST); // Carregando o problema
	
	for (size_t idx = 0; idx < node->forbidden_arcs.size(); idx++)
		p->cost[(*node).forbidden_arcs[idx].first][(*node).forbidden_arcs[idx].second] = 99999999;
	
	node->lower_bound = hungarian_solve(p);

	// cout << "Obj. value: " << obj_value << endl;
	
	// cout << "Assignment" << endl;
	// hungarian_print_assignment(p);
}

int main(int argc, char** argv) {

	Data* data = new Data(argc, argv[1]);
	data->readData();

	double **cost = new double*[data->getDimension()];
	for (int i = 0; i < data->getDimension(); i++) {
		cost[i] = new double[data->getDimension()];
		for (int j = 0; j < data->getDimension(); j++) {
			cost[i][j] = data->getDistance(i,j);
		}
	}

	Node root;
	list<Node> tree;
	tree.push_back(root);
	double upper_bound = 99999998; //numeric_limits<double>::infinity();
	double lower_bound = 0;
	int count = 0;
	int itmax = 100000;
	while (!tree.empty() && count < itmax) {
		count++;
		Node node = tree.front();
		tree.pop_front();
		hungarian_problem_t p;
		hungarian(&p, data, cost, &(node));
		findSubtours(&p, data, &(node));

		if (node.lower_bound > upper_bound) {
			// tree.erase(node);
			continue;
		}

		if (node.feasible) {
			upper_bound = min(upper_bound, node.lower_bound);
			// lower_bound = max(lower_bound, node.lower_bound);
			cerr << node.feasible << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";
		}
		else {
			for (size_t i = 0; i < node.subtours[node.chosen].size() - 1; i++) {
				Node n;
				n.forbidden_arcs = node.forbidden_arcs;
				std::pair<int,int> forbidden_arc = {
					node.subtours[node.chosen][i],
					node.subtours[node.chosen][i + 1]
				};
				n.forbidden_arcs.push_back(forbidden_arc);
				tree.push_back(n);
			}
			Node n;
			n.forbidden_arcs = node.forbidden_arcs;
			std::pair<int,int> forbidden_arc = {
				node.subtours[node.chosen][node.subtours[node.chosen].size() - 1],
				node.subtours[node.chosen][0]
			};
			n.forbidden_arcs.push_back(forbidden_arc);
			tree.push_back(n);
		}
		if (count % (itmax/100) == 0)
			cerr << node.feasible << " " << node.subtours[node.chosen].size() << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";
		// tree.erase(node);
	}
	cerr << count << endl;
	
	for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	delete [] cost;
	delete data;
	return 0;
}
