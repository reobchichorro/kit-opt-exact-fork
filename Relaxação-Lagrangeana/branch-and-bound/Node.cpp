#include "Node.h"

void Node::Solve(const vector<vector<double>>& og_cost, vector<vector<double>>& cost) {
	int count = 0;
	
	double epsilon = INITEPS;
	size_t k = 0;
	double mi = 0;
	vector<int> g(cost.size(), 2);
	double modG = 0;
	bool skip = false;
	bool a = false;
	bool b = false;
	bool c = false;

	lambda = curr_lambda;

	while (epsilon >= MINEPS && !feasible && !skip) {
		a = false;
		b = false;
		c = false;
		count++;
		feasible = true;
		for (size_t i = 0; i < cost.size(); i++) {
			for (size_t j = i+1; j < cost.size(); j++) {
				cost[i][j] = og_cost[i][j] - curr_lambda[i] - curr_lambda[j];
			}
		}

		for (auto it = forbidden_edges.cbegin(); it != forbidden_edges.cend(); it++) {
			cost[(*it).first][(*it).second] = MAXCOST;
		}
		Kruskal kr(cost, 1);
		double mst_sol = kr.MST(cost.size());
		double real_cost = 0;

		g = vector<int>(cost.size(), 2);
		for (size_t i = 0; i < kr.edges.size(); i++) {
			g[kr.edges[i].first]--;
			g[kr.edges[i].second]--;
			real_cost += og_cost[kr.edges[i].first][kr.edges[i].second];
			// cout << kr.edges[i].first << "," << kr.edges[i].second << " ";
		}
		// cout << "\n";
		
		modG = sqrt(accumulate(g.cbegin(), g.cend(), 0, [](int a, int b){ return a + b*b; }));
		feasible = modG < EPS;

		if (mst_sol > LB + EPS) {
			LB = mst_sol;
			k = 0;
			a = true;
			lambda = curr_lambda;
			edges = kr.edges;
			b = true;
			// if (curr_LB > LB) {
			// 	LB = curr_LB;
			// }
		}
		else {
			k++;
			if (k >= KMAX) {
				k = 0;
				epsilon /= 2;
			}
		}
		
		if (feasible) {
			if (real_cost + EPS < UB) {
				c = true;
				UB = min(UB, real_cost);
				LB = mst_sol;
				edges = kr.edges;
				lambda = curr_lambda;
				// lower_bound = min(lower_bound, node.lower_bound);
				// cerr << count << " c " << node.feasible << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";
			}
		}
		else {
			mi = epsilon*(UB - mst_sol)/modG;
	
			// cerr << mst_sol << "\t" << this_lower_bound << "\t" << lower_bound << "\t" << upper_bound << "\n";
			for (size_t i = 0; i < curr_lambda.size(); i++) {
				curr_lambda[i] += mi*g[i];
				if (curr_lambda[i] > MAXCOST/2 || curr_lambda[i] < -MAXCOST/2)
					skip = true;
				// cout << lambda[i] << " ";
			}
			// cout << "\n";
		}

		if (!b && c) {
			cerr << "aqui\n";
		}
	}
	// cerr << "haha\n";
	
	g = vector<int>(cost.size(), 0);
	for (size_t i = 0; i < edges.size(); i++) {
		g[edges[i].first]++;
		g[edges[i].second]++;
		real_cost += og_cost[edges[i].first][edges[i].second];
		// cout << edges[i].first << "," << edges[i].second << " ";
	}
	// cout << "\n";
	for (size_t i = 1; i < g.size(); i++) {
		if (g[i] > g[biggest])
			biggest = i;
	}

	// cout << LB << " " << UB << " " << count << " " << biggest << " " << g[biggest] << " " << forbidden_edges.size() << "\n";
}
