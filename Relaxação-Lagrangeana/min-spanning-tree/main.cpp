#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <numeric>

#include "data.h"
#include "Kruskal.h"

using namespace std;
#define EPS 1e-6
#define MAX_TREE_SIZE 1e+6

struct Node {
	double lower_bound; // cost of MST
	bool feasible = true;
};

double createInitialSolution(Data *data, const vector<vector<double>>& cost) {
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

int main(int argc, char** argv) {

	Data* data = new Data(argc, argv[1]);
	data->readData();

	// double **cost = new double*[data->getDimension()];
	// for (int i = 0; i < data->getDimension(); i++) {
	// 	cost[i] = new double[data->getDimension()];
	// 	for (int j = 0; j < data->getDimension(); j++) {
	// 		cost[i][j] = data->getDistance(i,j);
	// 		cerr << cost[i][j] << "\t";
	// 	}
	// 	cerr << "\n";
	// }
	
	vector<vector<double>> og_cost(data->getDimension(), vector<double>(data->getDimension()));
	vector<vector<double>> cost(data->getDimension(), vector<double>(data->getDimension()));
	for (int i = 0; i < data->getDimension(); i++) {
		for (int j = 0; j < data->getDimension(); j++) {
			og_cost[i][j] = data->getDistance(i,j);
			cost[i][j] = data->getDistance(i,j);
		}
	}
	
	// Node root;
	// list<Node> tree;
	// tree.push_back(root);
	double upper_bound = createInitialSolution(data, cost);//99999998; //numeric_limits<double>::infinity();
	cout << upper_bound << endl;
	double lower_bound = 0;

	int count = 0;
	long long int itmax = 100000000;
	// size_t curr_biggest = 0;
	// size_t max_forbidden = data->getDimension()*(data->getDimension()-1)/2 - data->getDimension();
	size_t tree_size = 1;

	vector<double> lambda(data->getDimension(), 0);
	double epsilon = 1;
	size_t k = 0;
	double mi = 0;
	vector<int> g(data->getDimension(), 2);
	double modG = 0;
	bool feasible = false;

	while (/*!tree.empty() && */count < itmax && tree_size < MAX_TREE_SIZE && (epsilon >= EPS && !feasible)) {
		// Node node = tree.back();
		// tree.pop_back();
		// tree_size--;
		// Sol solution;
		feasible = true;
		for (int i = 0; i < data->getDimension(); i++) {
			for (int j = 0; j < data->getDimension(); j++) {
				cost[i][j] = og_cost[i][j] - lambda[i] - lambda[j];
			}
		}

		Kruskal kr(cost, 1);
		double mst_sol = kr.MST(data->getDimension());
		double this_lower_bound = 0;

		g = vector<int>(data->getDimension(), 2);
		for (size_t i = 0; i < kr.edges.size(); i++) {
			g[kr.edges[i].first]--;
			g[kr.edges[i].second]--;
			this_lower_bound += og_cost[kr.edges[i].first][kr.edges[i].second];
		}

		modG = 0;
		for (size_t i = 0; i < g.size(); i++) {
			modG += g[i]*g[i];
		}
		double modGG = sqrt(accumulate(g.cbegin(), g.cend(), 0, [](int a, int b){ return a + b*b; }));

		mi = epsilon*(upper_bound - this_lower_bound)/modGG;

		for (size_t i = 0; i < lambda.size(); i++) {
			lambda[i] += mi*g[i];
			cerr << lambda[i] << " ";
		}
		int sum = accumulate(g.cbegin(), g.cend(), 0, [](int a, int b){ return a + abs(b); });
		cerr << sum << "\n";
		feasible = sum == 0;

		// if (count % 500000 == 0) {
		// 	// for (int i = 0; i < data->getDimension(); i++) {
		// 	// 	cerr << p.successors[i] << "\t";
		// 	// }
		// 	// cerr << "\n";
		// 	cerr << count << " c " << node.feasible << " " << node.chosen_subtour.size() << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";		// tree.erase(node);
		// }
		// count++;

		// if (this_lower_bound > upper_bound - EPS) {
		// 	// tree.erase(node);
		// 	break;
		// }

		if (mst_sol > lower_bound) {
			lower_bound = mst_sol;
			k = 0;
		}
		else {
			k++;
			if (k >= 10) {
				k = 0;
				epsilon /= 2;
			}
		}

		if (feasible) {
			if (this_lower_bound < upper_bound + EPS) {
				// for (int i = 0; i < data->getDimension(); i++) {
				// 	cerr << p.successors[i] << "\t";
				// }
				cerr << "\n\t";
				upper_bound = min(upper_bound, this_lower_bound);
				// lower_bound = min(lower_bound, node.lower_bound);
				// cerr << count << " c " << node.feasible << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";
			}
		}
		// else {
		// 	for (size_t i = 0; i < node.chosen_subtour.size() - 1; i++) {
		// 		Node n;

		// 		tree.push_back(n);
		// 	}
		// 	tree_size += node.chosen_subtour.size();// - 1;
		// 	// if (node.chosen_subtour.size() > 2 || true) {
		// 	Node n;

		// 	tree.push_back(n);
		// 	//tree_size++;
		// 	// }
		// }
	}
	cout << count << " " << upper_bound << " " << tree_size << endl;

	// for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	// delete [] cost;
	delete data;
	return 0;
}
