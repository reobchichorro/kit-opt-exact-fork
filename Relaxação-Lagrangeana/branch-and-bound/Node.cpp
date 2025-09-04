#include "Node.h"

Node::Node(long double lb, long double ub, const vector<long double>& _lambda, vector<long double>& _curr_lambda, vector<pair<size_t,size_t>>& _best_edges, vector<pair<size_t,size_t>>& _curr_edges) {
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

void Node::Solve(const vector<vector<long double>>& og_cost, vector<vector<long double>>& cost, vector <size_t>& _pset, vector<ii>& _non_zero_edges, vector<int>& g) {
	int count = 0;
	
	long double epsilon = INITEPS;
	size_t k = 0;
	long double mi = 0;
	long long int len2G = 0;
	bool keep_going = true;
	long double mst_sol = 0;

	for (size_t i=0; i< lambda.size(); i++)
		(*curr_lambda)[i] = lambda[i];

	// cout << "\tCnt\tG2\tMST\tRealCost\teps\tk\n";
	while (keep_going) {
		count++;
		// cerr << count << " count\n";
		feasible = true;
		mst_sol = 0;
		for (size_t i = 0; i < cost.size(); i++) {
			for (size_t j = i+1; j < cost.size(); j++) {
				cost[i][j] = og_cost[i][j] - (*curr_lambda)[i] - (*curr_lambda)[j];
			}
			g[i] = 2;
			mst_sol += 2*(*curr_lambda)[i];
		}

		for (auto it = forbidden_edges.cbegin(); it != forbidden_edges.cend(); it++) {
			cost[(*it).first][(*it).second] = MAXCOST;
		}
		Kruskal kr(cost, 1, _pset, _non_zero_edges, curr_edges);
		mst_sol += kr.MST(cost.size());

		long double real_cost = 0;

		for (size_t i = 0; i < (*curr_edges).size(); i++) {
			g[(*curr_edges)[i].first]--;
			g[(*curr_edges)[i].second]--;
			real_cost += og_cost[(*curr_edges)[i].first][(*curr_edges)[i].second];
			// if (count==1)
			// cerr << (*curr_edges)[i].first << "\t" << (*curr_edges)[i].second << "\n";
		}
		
		len2G = accumulate(g.cbegin(), g.cend(), 0, [](int a, int b){ return a + b*b; });
		feasible = len2G == 0;

		// if (len2G <= 2 && mst_sol > 674 && mst_sol < 676) {
		// 	cout << "\t" << count << "\t" << len2G << "\t" << mst_sol << "\t" << real_cost << "\t" << epsilon << "\t" << k << "\n";
		// 	cout << "\tLambda:\n\t\t";
		// 	for (size_t i=0; i< lambda.size(); i++)
		// 		cout << (*curr_lambda)[i] << "\t";
		// 	cout << "\n\t\t";
		// 	for (size_t i=0; i< lambda.size(); i++)
		// 		cout << lambda[i] << "\t";
		// 	cout << "\n";
		// }

		if (mst_sol >= UB) { // Might be wrong for TSP, since OBJ oscilates
			keep_going = false;
		}
		else if (mst_sol > LB) {// + EPS) {
			LB = max(mst_sol, LB);
			improved = LB < UB;
			k = 0;
			for (size_t i=0; i< lambda.size(); i++)
				lambda[i] = (*curr_lambda)[i];
			for (size_t i=0; i< (*best_edges).size(); i++)
				(*best_edges)[i] = (*curr_edges)[i];
		}
		else {
			k++;
			if (feasible)
				cerr << "aqui\n";
			if (k >= KMAX) {
				k = 0;
				epsilon /= 2;
			}
		}

		if (feasible) {
			if (abs(real_cost-mst_sol) > EPS) {
				cerr << count << "\t" << mst_sol << "\t" << real_cost << "\t" << LB << "\t" << UB << "\n";
			}
			if (real_cost + EPS < UB) {
				UB = min(UB, real_cost);
				improved = true;
				LB = max(mst_sol, LB);
				for (size_t i=0; i< lambda.size(); i++)
					lambda[i] = (*curr_lambda)[i];
				for (size_t i=0; i< (*best_edges).size(); i++)
					(*best_edges)[i] = (*curr_edges)[i];
			}
		}
		else {
			mi = epsilon*(UB - mst_sol)/(len2G);
			for (size_t i = 0; i < (*curr_lambda).size(); i++) {
				(*curr_lambda)[i] += mi*g[i];
				// if ((*curr_lambda)[i] < 0)
				// 	(*curr_lambda)[i] = 0;
				if ((*curr_lambda)[i] > MAXCOST/2 || (*curr_lambda)[i] < -MAXCOST/2)
					keep_going = false;
				// cerr << (*curr_lambda)[i] << "\n";
			}
		}

		if (epsilon < MINEPS || feasible)// || mst_sol >= UB)
			keep_going = false;
		// keep_going = epsilon >= MINEPS && !feasible && !skip && LB < UB;
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
		// cout << (*best_edges)[i].first << "," << (*best_edges)[i].second << " ";
	}
	// cout << "\n";
	// for (size_t i = 0; i < forbidden_edges.size(); i++) {
	// 	cout << forbidden_edges[i].first << "," << forbidden_edges[i].second << " ";
	// }
	// cout << "\n";
	for (size_t i = 1; i < g.size(); i++) {
		if (g[i] > g[biggest])
			biggest = i;
	}

	// cout << LB << " " << UB << " " << count << " " << biggest << " " << g[biggest] << " " << forbidden_edges.size() << "\n";
}
