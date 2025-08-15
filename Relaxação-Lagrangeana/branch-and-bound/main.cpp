#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <queue>

#include "data.h"
#include "Node.h"
#include "Kruskal.h"

using namespace std;
#define EPS 1e-6
#define MAX_TREE_SIZE 1e+8
#define BRANCHING 0 // 0 for DFS, 1 for BFS, 2 for LB
#define MAXCOST 99999999

double createInitialSolution(Data *data, const vector<vector<double>>& cost, vector<pair<size_t,size_t>>& edges) {
	unordered_set<int> remaining;
	list<size_t> sol; sol.push_back(0);
	for (int i = 1; i < data->getDimension(); i++)
		remaining.insert(i);
	
	double sol_cost = 0;
	size_t nearest = -1;
	double dist_nearest = 99999999;
	pair<size_t,size_t> curr_edge;
	while (!remaining.empty()) {
		for (auto it = remaining.begin(); it != remaining.end(); it++) {
			if (cost[sol.back()][*it] < dist_nearest) {
				nearest = *it;
				dist_nearest = cost[sol.back()][*it];
			}
		}
		edges.push_back(make_pair(sol.back(), nearest));
		sol.push_back(nearest);
		sol_cost += dist_nearest;
		remaining.erase(nearest);
		nearest = -1;
		dist_nearest = 99999999;
	}
	sol_cost += cost[sol.back()][0];
	for (size_t i = 0; i < edges.size(); i++) {
		cout << edges[i].first << "," << edges[i].second << " ";
	}
	cout << "\n";
	for (auto it = sol.cbegin(); it != sol.cend(); it++) {
		cout << *it << " ";
	}
	cout << sol_cost << " init\n";
	return sol_cost;
}

int main(int argc, char** argv) {

	Data* data = new Data(argc, argv[1]);
	data->readData();
	
	vector<vector<double>> og_cost(data->getDimension(), vector<double>(data->getDimension()));
	vector<vector<double>> cost(data->getDimension(), vector<double>(data->getDimension()));
	for (int i = 0; i < data->getDimension(); i++) {
		for (int j = 0; j < data->getDimension(); j++) {
			og_cost[i][j] = data->getDistance(i,j);
			cost[i][j] = data->getDistance(i,j);
		}
	}
	
	vector<pair<size_t,size_t>> edges;
	double UB = createInitialSolution(data, cost, edges);
	// double lower_bound = 0;
	cout << UB << endl;

	auto lambda = vector<double>(data->getDimension(), 0);
	Node root = Node(0, UB, lambda);

	list<Node> tree;
	tree.push_back(root);
	auto cmp = [](Node& a, Node&b){ return a.LB > b.LB; };
	priority_queue<Node, deque<Node>, decltype(cmp)> pq_tree(cmp);
	pq_tree.push(root);
	
	long long int count = 0;
	long long int itmax = 10000000;
	size_t tree_size = 1;

	size_t skipcount = 0;

	bool first = true;

	long long int count_lim = 5000;

	while ((BRANCHING == 2 ? !pq_tree.empty() : !tree.empty()) && count < itmax && tree_size < MAX_TREE_SIZE) {
		count++;
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

		if (UB - node.LB < 0.1) {
			skipcount++;
			continue;
		}
		
		node.Solve(og_cost, cost);

		if (count % count_lim == 0) {
			cerr << count << " " << skipcount << " c " << node.feasible << " " << node.LB << " " << UB << " " << node.forbidden_edges.size() << " " << tree_size << "\n";
			if (count > 10*count_lim)
				count_lim *= 20;
		}

		if (node.LB > UB || node.edges.size() == 0) {
			skipcount++;
			continue;
		}

		if (node.feasible) {
			if (node.UB < UB) {
				UB = min(UB, node.UB);
				// lower_bound = min(lower_bound, node.lower_bound);
				cerr << count << " " << skipcount << " c " << node.feasible << " " << node.LB << " " << UB << " " << node.forbidden_edges.size() << " " << tree_size << "\n";
			}
			if (first) {
				first = false;
				for (size_t i = 0; i < node.edges.size(); i++) {
					cerr << node.edges[i].first << "," << node.edges[i].second << " ";
				}
				cerr << "\n";
			}
		}

		// cerr << mst_sol << "\t" << this_lower_bound << "\t" << lower_bound << "\t" << UB << "\n";
		if (UB - node.LB < 0.1) {
			skipcount++;
			continue;
		}
		// cerr << "\n";

		// tree_size -= g[node.biggest]-2;
		for (size_t i = 0; i < node.edges.size(); i++) {
			if (node.biggest == node.edges[i].first || node.biggest == node.edges[i].second) {
				Node n(node.LB, UB, node.lambda);
				
				n.forbidden_edges = node.forbidden_edges;
				n.forbidden_edges.push_back({
					node.edges[i].first,
					node.edges[i].second
				});
				
				if (BRANCHING == 0 || BRANCHING == 1) {
					tree.push_back(n);
				}
				else if (BRANCHING == 2) {
					pq_tree.push(n);
				}
				tree_size++;
			}
		}
	}
	// cerr << "\n";
	cout << data->getInstanceName() << ";" << UB << " " << tree_size << " " << count << endl;

	// for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	// delete [] cost;
	delete data;
	return 0;
}
