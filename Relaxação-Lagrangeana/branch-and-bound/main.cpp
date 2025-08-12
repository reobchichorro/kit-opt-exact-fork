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
#define MINEPS 1e-6
#define MAX_TREE_SIZE 1e+8
#define KMAX 12
#define BRANCHING 2 // 0 for DFS, 1 for BFS, 2 for LB
#define MAXCOST 99999999

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
	
	Node root;
	root.lambda = vector<double>(data->getDimension(), 0);

	list<Node> tree;
	tree.push_back(root);
	auto cmp = [](Node& a, Node&b){ return a.lower_bound > b.lower_bound; };
	priority_queue<Node, deque<Node>, decltype(cmp)> pq_tree(cmp);
	pq_tree.push(root);
	
	double upper_bound = createInitialSolution(data, cost);
	cout << upper_bound << endl;
	double lower_bound = 0;

	long long int count = 0;
	long long int itmax = 10000000;
	size_t tree_size = 1;

	double mi = 0;
	vector<int> g(data->getDimension(), 2);
	double modG = 0;
	bool feasible = false;

	double node_cost;

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
		feasible = true;

		for (int i = 0; i < data->getDimension(); i++) {
			for (int j = 0; j < data->getDimension(); j++) {
				cost[i][j] = og_cost[i][j] - node.lambda[i] - node.lambda[j];
			}
			// cerr << node.lambda[i] << " ";
		}
		// cerr << "\n";

		for (auto it = node.forbidden_edges.cbegin(); it != node.forbidden_edges.cend(); it++) {
			cost[(*it).first][(*it).second] = MAXCOST;
			// cerr << (*it).first << "," << (*it).second << " ";
		}
		// cerr << "fa\n";

		Kruskal kr(cost, 1);
		node_cost = kr.MST(data->getDimension());

		g = vector<int>(data->getDimension(), 2);
		for (size_t i = 0; i < kr.edges.size(); i++) {
			g[kr.edges[i].first]--;
			g[kr.edges[i].second]--;
			node.real_cost += og_cost[kr.edges[i].first][kr.edges[i].second];
			// cerr << kr.edges[i].first << "," << kr.edges[i].second << " ";
		}
		// cerr << "edges\n" << node_cost << " " << node.real_cost << " " << node.lower_bound << " " << upper_bound << "\n";
		
		modG = sqrt(accumulate(g.cbegin(), g.cend(), 0, [](int a, int b){ return a + b*b; }));
		feasible = modG < EPS;

		if (count % 500000 == 0) {
			cerr << count << " c " << feasible << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_edges.size() << "\n";
		}

		if (node_cost > upper_bound || node.epsilon < MINEPS)
			continue;

		if (node_cost > node.lower_bound + EPS) {
			node.lower_bound = node_cost;
			node.k = 0;
		}
		else {
			node.k++;
			if (node.k >= KMAX) {
				node.k = 0;
				node.epsilon /= 2;
				// cerr << node.epsilon << "\n";
			}
		}

		if (feasible) {
			if (node.real_cost < upper_bound + EPS) {
				upper_bound = min(upper_bound, node.real_cost);
				// lower_bound = min(lower_bound, node.lower_bound);
				cerr << count << " c " << feasible << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_edges.size() << "\n";
			}
		}
		else {
			mi = node.epsilon*(upper_bound - node.lower_bound)/modG; // TODO: node.lower_bound vs node_cost
	
			// cerr << mst_sol << "\t" << this_lower_bound << "\t" << lower_bound << "\t" << upper_bound << "\n";
			bool skip = false;
			for (size_t i = 0; i < node.lambda.size(); i++) {
				node.lambda[i] += mi*g[i];
				if (node.lambda[i] > MAXCOST/10 || node.lambda[i] < -MAXCOST/10)
					skip = true;
				// cerr << node.lambda[i] << " ";
			}
			if (skip)
				continue;
			// cerr << "\n";

			for (size_t i = 1; i < g.size(); i++) {
				if (g[i] < g[node.biggest])
					node.biggest = i;
			}
			tree_size -= g[node.biggest]-2;
			for (size_t i = 0; i < kr.edges.size(); i++) {
				if (node.biggest == kr.edges[i].first || node.biggest == kr.edges[i].second) {
					Node n;
					n.lower_bound = node.lower_bound;
					
					n.forbidden_edges = node.forbidden_edges;
					n.forbidden_edges.push_back({
						kr.edges[i].first,
						kr.edges[i].second
					});
					
					n.lambda = node.lambda;
					n.epsilon = node.epsilon;
					n.k = node.k;
					
					if (BRANCHING == 0 || BRANCHING == 1) {
						tree.push_back(n);
					}
					else if (BRANCHING == 2) {
						pq_tree.push(n);
					}
				}
			}
			tree_size++;
			{
				Node n;
				n.lower_bound = node.lower_bound;
				
				n.forbidden_edges = node.forbidden_edges;
				n.lambda = node.lambda;
				n.epsilon = node.epsilon;
				n.k = node.k;
				
				if (BRANCHING == 0 || BRANCHING == 1) {
					tree.push_back(n);
				}
				else if (BRANCHING == 2) {
					pq_tree.push(n);
				}
			}
		}
	}
	// cerr << "\n";
	cout << data->getInstanceName() << ";" << upper_bound << " " << tree_size << " " << count << endl;

	// for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	// delete [] cost;
	delete data;
	return 0;
}
