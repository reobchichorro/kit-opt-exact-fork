#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>

#include "data.h"
#include "hungarian.h"
using namespace std;

#define EPS 1e-6
#define MAX_TREE_SIZE 1e+6

int n = -1;

struct Node {
	// vector<pair<int, int>> forbidden_arcs = vector<pair<int,int>>();
	vector<vector<bool>> forbidden_arcs;
	vector<vector<int>> subtours = vector<vector<int>>();
	double lower_bound; // cost of hungarian solution
	int chosen = -1; // index of the smallest subtours
	size_t biggest = 0; // size of the biggest subtours
	bool feasible = true;
};

// struct Sol {
// 	vector<size_t> sol = vector<size_t>(n);
// 	vector<vector<bool>> forbidden_arcs = vector<vector<bool>>(n, vector<bool>(n, 0));
// };

double createInitialSolution(Data *data, double **cost) {
	unordered_set<size_t> remaining;
	list<size_t> sol; sol.push_back(0);
	for (size_t i = 1; i < data->getDimension(); i++)
		remaining.insert(i);
	
	double sol_cost = 0;
	size_t nearest = -1;
	double dist_nearest = 99999999;
	while (!remaining.empty()) {
		for (auto it = remaining.begin(); it != remaining.end(); it++) {
			if (cost[sol.back()][*it] < dist_nearest) {
				nearest = *it;
				dist_nearest = cost[sol.back()][*it];
			}
		}
		sol.push_back(nearest);
		sol_cost += dist_nearest;
		remaining.erase(nearest);
		nearest = -1;
		dist_nearest = 99999999;
	}
	sol_cost += cost[sol.back()][0];
	for (auto it = sol.cbegin(); it != sol.cend(); it++) {
		cerr << *it << " ";
	}
	cerr << sol_cost << " init\n";
	return sol_cost;
}

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
		node->biggest = max((*it).second, node->biggest);
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
	
	// for (size_t idx = 0; idx < node->forbidden_arcs.size(); idx++) {
	// 	p->cost[(*node).forbidden_arcs[idx].first][(*node).forbidden_arcs[idx].second] = 99999999;
	// 	p->cost[(*node).forbidden_arcs[idx].second][(*node).forbidden_arcs[idx].first] = 99999999;
	// }

	for (int i = 0; i < data->getDimension(); i++) {
		for (int j = 0; j < data->getDimension(); j++) {
			p->cost[i][j] = node->forbidden_arcs[i][j]*99999999 + (1-node->forbidden_arcs[i][j])*p->cost[i][j];
		}
	}
	
	node->lower_bound = hungarian_solve(p);
	
	// solution.forbidden_arcs = node->forbidden_arcs;
	// solution.sol = vector<size_t>(data->getDimension());
	// for (int i = 0; i < data->getDimension(); i++)
	// 	solution.sol[i] = p->successors[i];
	// cout << "Obj. value: " << obj_value << endl;
	
	// cout << "Assignment" << endl;
	// hungarian_print_assignment(p);
}

int main(int argc, char** argv) {

	Data* data = new Data(argc, argv[1]);
	data->readData();
	n = data->getDimension();

	double **cost = new double*[data->getDimension()];
	for (int i = 0; i < data->getDimension(); i++) {
		cost[i] = new double[data->getDimension()];
		for (int j = 0; j < data->getDimension(); j++) {
			cost[i][j] = data->getDistance(i,j);
			cerr << cost[i][j] << "\t";
		}
		cerr << "\n";
	}

	Node root;
	root.forbidden_arcs = vector<vector<bool>>(data->getDimension(), vector<bool>(data->getDimension(), 0));
	list<Node> tree;
	tree.push_back(root);
	double upper_bound = createInitialSolution(data, cost);//99999998; //numeric_limits<double>::infinity();
	double lower_bound = root.lower_bound;

	// unordered_set<Sol> visited;

	int count = 0;
	long long int itmax = 10000000;
	size_t curr_biggest = 0;
	size_t max_forbidden = data->getDimension()*(data->getDimension()-1)/2 - data->getDimension();
	size_t tree_size = 1;

	while (!tree.empty() && count < itmax && tree_size < MAX_TREE_SIZE) {
		count++;
		Node node = tree.back();
		tree.pop_back();
		tree_size--;
		hungarian_problem_t p;
		// Sol solution;
		hungarian(&p, data, cost, &(node));

		findSubtours(&p, data, &(node));

		if (count % 100000 == 0)
			cerr << count << " c " << node.feasible << " " << node.subtours[node.chosen].size() << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";		// tree.erase(node);

		if (node.lower_bound > upper_bound - EPS) {
			// tree.erase(node);
			continue;
		}

		// if (node.biggest > curr_biggest) {
		// 	curr_biggest = node.biggest;
		// 	cerr << node.feasible << " " << node.subtours[node.chosen].size() << " " << curr_biggest << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";
		// }
		if (node.forbidden_arcs.size() > max_forbidden)
			cerr << "Forbidden: " << node.feasible << " " << node.subtours[node.chosen].size() << " " << curr_biggest << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";

		if (node.feasible) {
			upper_bound = min(upper_bound, node.lower_bound);
			// lower_bound = min(lower_bound, node.lower_bound);
			cerr << node.feasible << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";
		}
		else {
			for (size_t i = 0; i < node.subtours[node.chosen].size() - 1; i++) {
				Node n;
				n.forbidden_arcs = node.forbidden_arcs;
				n.forbidden_arcs[
					node.subtours[node.chosen][i]
				][
					node.subtours[node.chosen][i + 1]
				] = 1;
				n.forbidden_arcs[
					node.subtours[node.chosen][i + 1]
				][
					node.subtours[node.chosen][i]
				] = 1;
				// n.forbidden_arcs.push_back({
				// 	node.subtours[node.chosen][i],
				// 	node.subtours[node.chosen][i + 1]
				// });
				// n.forbidden_arcs.push_back({
				// 	node.subtours[node.chosen][i + 1],
				// 	node.subtours[node.chosen][i]
				// });
				tree.push_back(n);
			}
			Node n;
			n.forbidden_arcs = node.forbidden_arcs;
			n.forbidden_arcs[
				node.subtours[node.chosen][node.subtours[node.chosen].size() - 1]
			][
				node.subtours[node.chosen][0]
			] = 1;
			n.forbidden_arcs[
				node.subtours[node.chosen][0]
			][
				node.subtours[node.chosen][node.subtours[node.chosen].size() - 1]
			] = 1;
			// n.forbidden_arcs.push_back({
			// 	node.subtours[node.chosen][node.subtours[node.chosen].size() - 1],
			// 	node.subtours[node.chosen][0]
			// });
			// n.forbidden_arcs.push_back({
			// 	node.subtours[node.chosen][0],
			// 	node.subtours[node.chosen][node.subtours[node.chosen].size() - 1]
			// });
			tree.push_back(n);
			tree_size += node.subtours[node.chosen].size();
		}
	}
	cerr << count << " " << upper_bound << " " << tree_size << endl;
	
	for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	delete [] cost;
	delete data;
	return 0;
}
