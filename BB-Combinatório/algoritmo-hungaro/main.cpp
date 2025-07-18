#include <iostream>
using namespace std;

#include "data.h"
#include "hungarian.h"
#include <vector>
#include <list>

struct Node {
	vector<pair<int, int>> forbidden_arcs;
	vector<vector<int>> subtour;
	double lower_bound; // cost of hungarian solution
	int chosen; // index of the smallest subtours
	bool feasible;
};

void findSubtours(hungarian_problem_t *p, Data *data, Node* node) {
	
}

void hungarian(hungarian_problem_t *p, Data *data, double **cost, Node* node) {
	hungarian_init(p, cost, data->getDimension(), data->getDimension(), HUNGARIAN_MODE_MINIMIZE_COST); // Carregando o problema
	
	for (size_t idx = 0; idx < node->forbidden_arcs.size(); idx++)
		p->cost[(*node).forbidden_arcs[idx].first][(*node).forbidden_arcs[idx].second] = 99999999;
	
	double obj_value = hungarian_solve(p);
	// cout << "Obj. value: " << obj_value << endl;
	
	// cout << "Assignment" << endl;
	// hungarian_print_assignment(p);

	hungarian_free(p);
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
	double upper_bound = numeric_limits<double>::infinity();

	while (!tree.empty()) {
		auto node = tree.begin();
		hungarian_problem_t p;
		vector<vector<int>> subtour = getSolutionHungarian(*node);

		if (node->lower_bound > upper_bound) {
			tree.erase(node);
			continue;
		}

		if (node->feasible)
			upper_bound = min(upper_bound, node->lower_bound);
		else {
			for (int i = 0; i < node.subtour[root.chosen].size() - 1; i++) {
				Node n;
				n.arcos_proibidos = raiz.arcos_proibidos;
				std::pair<int,int> forbidden_arc = {
					node.subtour[root.chosen][i],
					node.subtour[root.chosen][i + 1]
				};
				n.forbidden_arcs.push_back(forbidden_arc);
				tree.push_back(n);
			}
		}
	}
	
	for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	delete [] cost;
	delete data;
	return 0;
}
