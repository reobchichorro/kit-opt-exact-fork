#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <queue>

#include "data.h"
#include "hungarian.h"

using namespace std;
#define EPS 1e-6
#define MAX_TREE_SIZE 1e+6
#define BRANCHING 2 // 0 for DFS, 1 for BFS, 2 for LB

int n = -1;

struct Node {
	vector<pair<int, int>> forbidden_arcs = vector<pair<int,int>>();
	// vector<vector<bool>> forbidden_arcs;
	vector<int> chosen_subtour = vector<int>(); // smallest subtour
	double lower_bound; // cost of hungarian solution
	// int chosen = -1; // index of the smallest subtours
	size_t biggest = 0; // size of the biggest subtours
	bool feasible = true;
};

// struct Sol {
// 	vector<size_t> sol = vector<size_t>(n);
// 	vector<vector<bool>> forbidden_arcs = vector<vector<bool>>(n, vector<bool>(n, 0));
// };

double createInitialSolution(Data *data, double **cost) {
	unordered_set<int> remaining;
	list<size_t> sol; sol.push_back(0);
	for (int i = 1; i < data->getDimension(); i++)
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
	vector<size_t> disjointedSets = vector<size_t>(data->getDimension());
	size_t j;
	for (size_t i = 0; i < disjointedSets.size(); i++)
		disjointedSets[i] = i;
	for (size_t i = 0; i < disjointedSets.size(); i++) {
		if (disjointedSets[i] != i)
			continue;
		j = p->successors[i];
		while(j != i) {
			disjointedSets[j] = disjointedSets[i];
			j = p->successors[j];
		}
	}
	unordered_map<size_t, size_t> subtoursSizes;
	for (size_t i = 0; i < disjointedSets.size(); i++) {
		if (subtoursSizes.contains(disjointedSets[i]))
			subtoursSizes[disjointedSets[i]]++;
		else
			subtoursSizes.emplace(disjointedSets[i], 1);
	}
	int chosen = -1;
	size_t smallest = disjointedSets.size() + 1;
	for (auto it = subtoursSizes.cbegin(); it != subtoursSizes.cend(); it++) {
		if ((*it).second < smallest) {
			chosen = (*it).first;
			smallest = (*it).second;
		}
	}
	
	int idx = chosen;
	do {
		node->chosen_subtour.push_back(idx);
		idx = p->successors[idx];

	} while(idx != chosen);
	if (smallest < disjointedSets.size())
		node->feasible = false;
}

void hungarian(hungarian_problem_t *p, Data *data, double **cost, Node* node) {
	hungarian_reset(p, cost, data->getDimension(), data->getDimension(), HUNGARIAN_MODE_MINIMIZE_COST); // Carregando o problema
	
	for (size_t idx = 0; idx < node->forbidden_arcs.size(); idx++) {
		p->cost[(*node).forbidden_arcs[idx].first][(*node).forbidden_arcs[idx].second] = 99999999;
		// p->cost[(*node).forbidden_arcs[idx].second][(*node).forbidden_arcs[idx].first] = 99999999;
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
	root.forbidden_arcs = vector<pair<int,int>>(); //vector<vector<bool>>(data->getDimension(), vector<bool>(data->getDimension(), 0));

	list<Node> tree;
	tree.push_back(root);
	auto cmp = [](Node& a, Node&b){ return a.lower_bound > b.lower_bound; };
	priority_queue<Node, deque<Node>, decltype(cmp)> pq_tree(cmp);
	pq_tree.push(root);
	
	double upper_bound = createInitialSolution(data, cost);//99999998; //numeric_limits<double>::infinity();
	cout << upper_bound << endl;
	hungarian_problem_t p;
	hungarian_init(&p, cost, data->getDimension(), data->getDimension(), HUNGARIAN_MODE_MINIMIZE_COST); // Carregando o problema

	int count = 0;
	long long int itmax = 100000000;
	size_t tree_size = 1;

	while ((BRANCHING == 2 ? !pq_tree.empty() : !tree.empty()) && count < itmax && tree_size < MAX_TREE_SIZE) {
		Node node;
		if (BRANCHING == 0) {
			node = tree.back();
			tree.pop_back();
		}
		else if (BRANCHING == 1) {
			node = tree.front();
			tree.pop_front();
		}
		else if (BRANCHING == 2) {
			node = pq_tree.top();
			pq_tree.pop();
		}
		tree_size--;
		// Sol solution;
		hungarian(&p, data, cost, &(node));

		findSubtours(&p, data, &(node));

		if (count % 500000 == 0) {
			for (int i = 0; i < data->getDimension(); i++) {
				cerr << i << "," << p.successors[i] << "\t";
			}
			cerr << "\n";
			cerr << count << " c " << node.feasible << " " << node.chosen_subtour.size() << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";		// tree.erase(node);
		}
		count++;

		if (node.lower_bound > upper_bound - EPS) {
			// tree.erase(node);
			continue;
		}

		if (node.feasible) {
			if (node.lower_bound < upper_bound + EPS) {
				for (int i = 0; i < data->getDimension(); i++) {
					cerr << i << "," << p.successors[i] << "\t";
				}
				cerr << "\n\t";
				upper_bound = min(upper_bound, node.lower_bound);
				// lower_bound = min(lower_bound, node.lower_bound);
				cerr << count << " c " << node.feasible << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";
			}
		}
		else {
			for (size_t i = 0; i < node.chosen_subtour.size() - 1; i++) {
				Node n;
				n.lower_bound = node.lower_bound;
				n.forbidden_arcs = node.forbidden_arcs;
				n.forbidden_arcs.push_back({
					node.chosen_subtour[i],
					node.chosen_subtour[i + 1]
				});
				if (BRANCHING == 0 || BRANCHING == 1) {
					tree.push_back(n);
				}
				else if (BRANCHING == 2) {
					pq_tree.push(n);
				}
			}
			tree_size += node.chosen_subtour.size();
			Node n;
			n.lower_bound = node.lower_bound;
			n.forbidden_arcs = node.forbidden_arcs;
			n.forbidden_arcs.push_back({
				node.chosen_subtour[node.chosen_subtour.size() - 1],
				node.chosen_subtour[0]
			});
			if (BRANCHING == 0 || BRANCHING == 1) {
				tree.push_back(n);
			}
			else if (BRANCHING == 2) {
				pq_tree.push(n);
			}
		}
	}
	cout << count << " " << upper_bound << " " << tree_size << endl;
	hungarian_free(&p);
	
	for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	delete [] cost;
	delete data;
	return 0;
}
