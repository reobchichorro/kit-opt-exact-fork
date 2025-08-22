#include "Node.h"

Node::Node(double lb, double ub, const vector<double>& _lambda, vector<double>& _curr_lambda, vector<pair<size_t,size_t>>& _best_edges, vector<pair<size_t,size_t>>& _curr_edges) {
	LB = lb;
	UB = ub;
	lambda = _lambda;
	curr_lambda = &_curr_lambda;
	for (size_t i=0; i< _lambda.size(); i++)
		(*curr_lambda)[i] = _lambda[i];
	best_edges = &_best_edges;
	curr_edges = &_curr_edges;
}

Node::~Node(){
	curr_lambda = NULL;
	best_edges = NULL;
	curr_edges = NULL;
}

void Node::Solve(const vector<vector<double>>& og_cost, vector<vector<double>>& cost, vector <size_t>& _pset, vector<ii>& _non_zero_edges, vector<int>& g) {
	int count = 0;
	
	double epsilon = INITEPS;
	size_t k = 0;
	double mi = 0;
	double modG = 0;
	bool skip = false;

	for (size_t i=0; i< lambda.size(); i++)
		(*curr_lambda)[i] = lambda[i];

	while (epsilon >= MINEPS && !feasible && !skip) {
		count++;
		feasible = true;
		for (size_t i = 0; i < cost.size(); i++) {
			for (size_t j = i+1; j < cost.size(); j++) {
				cost[i][j] = og_cost[i][j] - (*curr_lambda)[i] - (*curr_lambda)[j];
			}
			g[i] = 2;
		}

		for (auto it = forbidden_edges.cbegin(); it != forbidden_edges.cend(); it++) {
			cost[(*it).first][(*it).second] = MAXCOST;
		}
		Kruskal kr(cost, 1, _pset, _non_zero_edges, curr_edges);
		double mst_sol = kr.MST(cost.size());
		double real_cost = 0;

		for (size_t i = 0; i < (*curr_edges).size(); i++) {
			g[(*curr_edges)[i].first]--;
			g[(*curr_edges)[i].second]--;
			real_cost += og_cost[(*curr_edges)[i].first][(*curr_edges)[i].second];
			// cout << (*edges)[i].first << "," << (*edges)[i].second << " ";
		}
		// cout << "\n";
		
		modG = sqrt(accumulate(g.cbegin(), g.cend(), 0, [](int a, int b){ return a + b*b; }));
		feasible = modG < EPS;

		if (mst_sol > LB + EPS) {
			improved = true;
			LB = mst_sol;
			k = 0;
			for (size_t i=0; i< lambda.size(); i++)
				lambda[i] = (*curr_lambda)[i];
			for (size_t i=0; i< (*best_edges).size(); i++)
				(*best_edges)[i] = (*curr_edges)[i];
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
				UB = min(UB, real_cost);
				improved = true;
				LB = mst_sol;
				for (size_t i=0; i< lambda.size(); i++)
					lambda[i] = (*curr_lambda)[i];
				for (size_t i=0; i< (*best_edges).size(); i++)
					(*best_edges)[i] = (*curr_edges)[i];
				// lower_bound = min(lower_bound, node.lower_bound);
				// cerr << count << " c " << node.feasible << " " << node.lower_bound << " " << upper_bound << " " << node.forbidden_arcs.size() << "\n";
			}
		}
		else {
			mi = epsilon*(UB - mst_sol)/modG;
	
			// cerr << mst_sol << "\t" << this_lower_bound << "\t" << lower_bound << "\t" << upper_bound << "\n";
			for (size_t i = 0; i < (*curr_lambda).size(); i++) {
				(*curr_lambda)[i] += mi*g[i];
				if ((*curr_lambda)[i] > MAXCOST/2 || (*curr_lambda)[i] < -MAXCOST/2)
					skip = true;
				// cout << lambda[i] << " ";
			}
			// cout << "\n";
		}
	}
	// cerr << "haha\n";
	if (!improved)
		return;
	
	for (size_t i = 0; i < cost.size(); i++) {
		g[i] = 0;
	}
	for (size_t i = 0; i < (*best_edges).size(); i++) {
		g[(*best_edges)[i].first]++;
		g[(*best_edges)[i].second]++;
		real_cost += og_cost[(*best_edges)[i].first][(*best_edges)[i].second];
		// cout << edges[i].first << "," << edges[i].second << " ";
	}
	// cout << "\n";
	for (size_t i = 1; i < g.size(); i++) {
		if (g[i] > g[biggest])
			biggest = i;
	}

	// cout << LB << " " << UB << " " << count << " " << biggest << " " << g[biggest] << " " << forbidden_edges.size() << "\n";
}
