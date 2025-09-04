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
#define BRANCHING 2 // 0 for DFS, 1 for BFS, 2 for LB
#define MAXCOST 99999999

auto cmp = [](Node& a, Node&b){ return a.LB > b.LB; };

void setUBs(unordered_map<string, long double>& UBs) {
	UBs["a280"] = 2579;
	UBs["ali535"] = 202384;
	UBs["att48"] = 10628;
	UBs["att532"] = 27731;
	UBs["bayg29"] = 1610;
	UBs["bays29"] = 2020;
	UBs["berlin52"] = 7542;
	UBs["bier127"] = 118282;
	UBs["brazil58"] = 25395;
	UBs["brg180"] = 1950;
	UBs["burma14"] = 3323;
	UBs["ch130"] = 6110;
	UBs["ch150"] = 6528;
	UBs["d198"] = 15780;
	UBs["d493"] = 35042;
	UBs["dantzig42"] = 699;
	UBs["eil101"] = 629;
	UBs["eil51"] = 426;
	UBs["eil76"] = 538;
	UBs["fl417"] = 11861;
	UBs["fri26"] = 937;
	UBs["gil262"] = 2378.7;
	UBs["gr120"] = 6942;
	UBs["gr137"] = 69853;
	UBs["gr17"] = 2085;
	UBs["gr202"] = 40160.1;
	UBs["gr21"] = 2707;
	UBs["gr229"] = 134613;
	UBs["gr24"] = 1272;
	UBs["gr431"] = 171530;
	UBs["gr48"] = 5046;
	UBs["gr96"] = 55209;
	UBs["hk48"] = 11461;
	UBs["kroA100"] = 21282;
	UBs["kroA150"] = 26524;
	UBs["kroA200"] = 29368;
	UBs["kroB100"] = 22141;
	UBs["kroB150"] = 26130;
	UBs["kroB200"] = 29437.2;
	UBs["kroC100"] = 20749;
	UBs["kroD100"] = 21294;
	UBs["kroE100"] = 22068;
	UBs["lin105"] = 14379;
	UBs["lin318"] = 42045.7;
	UBs["linhp318"] = 42053.1;
	UBs["pcb442"] = 50876;
	UBs["pr107"] = 44303;
	UBs["pr124"] = 59030;
	UBs["pr136"] = 96772;
	UBs["pr144"] = 58537;
	UBs["pr152"] = 73682;
	UBs["pr226"] = 80369;
	UBs["pr264"] = 49135;
	UBs["pr299"] = 48194.8;
	UBs["pr76"] = 108159;
	UBs["rat195"] = 2326.1;
	UBs["rat99"] = 1211;
	UBs["rd100"] = 7910;
	UBs["rd400"] = 15296.1;
	UBs["si175"] = 21407;
	UBs["si535"] = 48466.8;
	UBs["st70"] = 675;
	UBs["swiss42"] = 1273;
	UBs["ts225"] = 126643;
	UBs["tsp225"] = 3916;
	UBs["u159"] = 42080;
	UBs["ulysses16"] = 6859;
	UBs["ulysses22"] = 7013;
	for (auto it = UBs.begin(); it != UBs.end(); it++)
		UBs[it->first] = it->second + 1;
}

long double createInitialSolution(Data *data, const vector<vector<long double>>& cost, vector<pair<size_t,size_t>>& edges) {
	unordered_set<int> remaining;
	list<size_t> sol; sol.push_back(0);
	for (int i = 1; i < data->getDimension(); i++)
		remaining.insert(i);
	
	long double sol_cost = 0;
	size_t nearest = -1;
	long double dist_nearest = 99999999;
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
	edges.push_back(make_pair(sol.back(), 0));
	sol_cost += cost[sol.back()][0];
	// for (size_t i = 0; i < edges.size(); i++) {
	// 	cout << edges[i].first << "," << edges[i].second << " ";
	// }
	// cout << "\n";
	// for (auto it = sol.cbegin(); it != sol.cend(); it++) {
	// 	cout << *it << " ";
	// }
	// cout << sol_cost << " init\n";
	return sol_cost;
}

long double getLB(const list<Node>& tree, const priority_queue<Node, deque<Node>, decltype(cmp)>& pq_tree, long double feas_lb) {
	long double LB = MAXCOST;
	if (BRANCHING == 2) {
		if (!pq_tree.empty())
			LB = pq_tree.top().LB;
	}
	else {
		for (auto it = tree.cbegin(); it != tree.cend(); it++)
			LB = min(LB, it->LB);
	}
	return min(LB, feas_lb);
}

int main(int argc, char** argv) {

	Data* data = new Data(argc, argv[1]);
	data->readData();

	unordered_map<string, long double> UBs;
	setUBs(UBs);
	
	vector <size_t> pset(data->getDimension());
	
	vector<vector<long double>> og_cost(data->getDimension(), vector<long double>(data->getDimension()));
	vector<vector<long double>> cost(data->getDimension(), vector<long double>(data->getDimension()));
	for (int i = 0; i < data->getDimension(); i++) {
		for (int j = 0; j < data->getDimension(); j++) {
			og_cost[i][j] = data->getDistance(i,j);
			cost[i][j] = data->getDistance(i,j);
			// cerr << og_cost[i][j] << "\t";
		}
		// cerr << "\n";
	}
	
	vector<int> g(cost.size());
	vector<ii> non_zero_edges;
	for (size_t i = 1; i < cost.size(); ++i) {
		for (size_t j = i+1; j < cost[i].size(); ++j) {
			non_zero_edges.push_back( make_pair(i, j) );
			// cerr << i << "\t" << j << "\n";
		}
	}
	
	vector<pair<size_t,size_t>> best_edges;
	vector<pair<size_t,size_t>> node_best_edges;
	vector<pair<size_t,size_t>> node_curr_edges;
	long double UB = min(createInitialSolution(data, cost, best_edges), UBs[data->getInstanceName()]);
	node_best_edges = node_curr_edges = best_edges;
	cout << "Initial UB: " << UB << endl;
	long double LB = 0; // Best LB
	long double feas_lb = MAXCOST; // LB of current best feasible solution

	auto lambda = vector<long double>(data->getDimension(), 0);
	auto curr_lambda = vector<long double>(data->getDimension(), 0);
	Node root = Node(0, UB, lambda, curr_lambda, node_best_edges, node_curr_edges);

	list<Node> tree;
	tree.push_back(root);	
	priority_queue<Node, deque<Node>, decltype(cmp)> pq_tree(cmp);
	pq_tree.push(root);
	
	long long int count = 0;
	long long int itmax = 100000;
	size_t tree_size = 1;

	size_t skipcount = 0;

	long long int count_lim = 250;
	bool updateLB = false;

	cout << "Visited Nodes\tSkipped Nodes\tLB\tUB\tTree Size\tGap\n";

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
		// bool stop = node.code.size() == 1 && node.code[0] != 'c';

		if (UB < node.LB + EPS) {
			skipcount++;
			// cout << count << "\t" << skipcount << "\t" << LB << "\t" << UB << "\t" << tree_size << "\t" << (UB-LB)/LB << "\t" << node.LB << "\tSkipped\n";
			continue;
		}
		updateLB = node.LB == LB;

		node.UB = UB;
		node.Solve(og_cost, cost, pset, non_zero_edges, g);

		if (count % count_lim == 0) {
			// LB = getLB(tree, pq_tree);
			cout << count << "\t" << skipcount << "\t" << LB << "\t" << UB << "\t" << tree_size << "\t" << (UB-LB)/LB << endl;
			if (count >= 10*count_lim)
				count_lim *= 2;
		}

		if (node.feasible) {
			if (node.UB + EPS < UB) {
				UB = min(UB, node.UB);
				feas_lb = node.LB;
				LB = getLB(tree, pq_tree, feas_lb);
				updateLB = false;
				cout << count << "\t" << skipcount << "\t" << LB << "\t" << UB << "\t" << tree_size << "\t" << (UB-LB)/LB << "\tImproved solution found" << "\n";
				// if (true) {
				// 	// first = false;
				// 	for (size_t i = 0; i < node_best_edges.size(); i++) {
				// 		cout << node_best_edges[i].first << "," << node_best_edges[i].second << " ";
				// 	}
				// 	cout << "\n";
				// 	for (size_t i = 0; i < node.forbidden_edges.size(); i++) {
				// 		cout << node.forbidden_edges[i].first << "," << node.forbidden_edges[i].second << " ";
				// 	}
				// 	cout << "\n";
				// }
				for (size_t i=0; i< best_edges.size(); i++)
					best_edges[i] = (*node.best_edges)[i];
			}
		}
		else if (UB < node.LB + EPS || !node.improved) {
			if (updateLB)
				LB = getLB(tree, pq_tree, feas_lb);
			skipcount++;
			// if (node.improved)
			// 	cout << count << "\t" << skipcount << "\t" << LB << "\t" << UB << "\t" << tree_size << "\t" << (UB-LB)/LB << "\t" << node.LB << "\tDead end\n";
			// else
			// 	cout << count << "\t" << skipcount << "\t" << LB << "\t" << UB << "\t" << tree_size << "\t" << (UB-LB)/LB << "\tNot improved\n";
			continue;
		} else {
			// char suffix = 'a';
			for (size_t i = 0; i < node_best_edges.size(); i++) {
				if (node.biggest == node_best_edges[i].first || node.biggest == node_best_edges[i].second) {
					Node n(node.LB, UB, node.lambda, curr_lambda, node_best_edges, node_curr_edges);
	
					n.forbidden_edges = node.forbidden_edges;
					n.forbidden_edges.push_back({
						node_best_edges[i].first,
						node_best_edges[i].second
					});
					
					if (BRANCHING == 0 || BRANCHING == 1) {
						tree.push_back(n);
					}
					else if (BRANCHING == 2) {
						pq_tree.push(n);
					}
					tree_size++;
					// suffix++;
				}
			}
		}

		if (updateLB)
			LB = getLB(tree,pq_tree, feas_lb);
	}
	LB = getLB(tree, pq_tree, feas_lb);
	
	if (tree_size != 0) {
		cout << "\nBest solution found for " << data->getInstanceName() << "\tUB\tLB\tGap\tVisited Nodes\tTree Size\n";
		cout << "\t" << UB << "\t" << LB << "\t" << (UB-LB)/LB << count << "\t" << tree_size << "\n";
	}
	else if (UB - LB < EPS) {
		cout << "\nOptimal found for " << data->getInstanceName() << "\tOPT\tVisited Nodes\n";
		cout << "\t" << UB << "\t" << count << "\n";
	} else {
		cout << "\nOptimal found? for " << data->getInstanceName() << "\tOPT\tLB\tVisited Nodes\n";
		cout << "\t" << UB << "\t" << LB << "\t" << count << "\n";
	}


	cout << "Best solution edges:\n";
	for (size_t i = 0; i < node_best_edges.size(); i++) {
		cout << best_edges[i].first << "," << best_edges[i].second << " ";
	}
	cout << "\n";

	// for (int i = 0; i < data->getDimension(); i++) delete [] cost[i];
	// delete [] cost;
	delete data;
	return 0;
}
