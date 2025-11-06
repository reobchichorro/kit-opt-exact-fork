#include <ilcplex/ilocplex.h>
#include "combo.c"
#include <vector>
#include <cstdlib>
#include <list>
#include <unordered_set>

using namespace std;

#define EPS 1e-6

bool AllInt(const IloCplex& rmp, const IloNumVarArray& lambda, const vector<vector<bool>>& columns, const int n, vector<vector<double>>& pairs, size_t& fraci, size_t& fracj) {
	bool all_int = true;
	vector<double> values(n, 0);

	for (size_t k = 0; k < lambda.getSize(); k++)
	{
		if (rmp.getValue(lambda[k]) > EPS && abs(1 - rmp.getValue(lambda[k])) > EPS) {
			all_int = false;
			break;
		}
	}

	if (all_int)
		return all_int;

	for (size_t i = 0; i < n; i++) {
		for (size_t j = i+1; j < n; j++) {
			pairs[i][j] = 0;
			for (size_t k = 0; k < lambda.getSize(); k++)
			{
				if (columns[k][i] && columns[k][j]) {
					pairs[i][j] += rmp.getValue(lambda[k]);
				}
			}
		}
	}

	double frac = 2;

	for (size_t i = 0; i < n; i++) {
		for (size_t j = i+1; j < n; j++) {
			// cout << pairs[i][j] << " ";
			if (abs(0.5 - pairs[i][j]) + EPS < abs(0.5 - frac)) {
				fraci = i; fracj = j;
				frac = pairs[i][j];
			}
		}
		// cout << "\n";
	}

	return all_int;
}

struct Node {
	unordered_set<size_t> used_columns;
	vector<pair<pair<size_t, size_t>, bool>> pairs;
	double LB;
};

void solveMIP() {

}

int main() 
{
	const double M = 1e7;
	int n = 5;
	vector<int> weight = {2, 1, 3, 3, 5};
	int capacity = 7;
	bool use_mip = false;

	cin >> n >> capacity;
	weight = vector<int>(n);
	for (int i=0; i < n; i++)
	cin >> weight[i];

	long long int UB = 9999999;

	vector<vector<bool>> columns(n, vector<bool>(n,0));

	item* items = new item[n];
	for (int i = 0; i < n; i++)
	{
		items[i].w = weight[i];
		items[i].p = M;
		items[i].x = 0;
		items[i].index = i;
		columns[i][i] = 1;
	}

	IloEnv env;
	IloModel master_model(env);

	IloNumVarArray lambda(env, n, 0, IloInfinity);

	IloExpr sum_obj(env);
	IloRangeArray partition_constraint(env);

	for (int i = 0; i < n; i++)
	{
		char var_name[50];
		sprintf(var_name, "y%d", i);

		lambda[i].setName(var_name);
		sum_obj += M * lambda[i];

		partition_constraint.add(lambda[i] == 1);
	}

	master_model.add(partition_constraint);

	IloObjective master_objective = IloMinimize(env, sum_obj);
	master_model.add(master_objective);

	IloCplex rmp(master_model);

	rmp.setOut(env.getNullStream()); // disables CPLEX log

	rmp.solve();

	cout << "Initial lower bound: " << rmp.getObjValue() << endl;

	cout << "Initial solution: ";
	for (size_t j = 0; j < lambda.getSize(); j++)
	{
		cout << rmp.getValue(lambda[j]) << " ";
	}
	cout << endl;

	int lambda_counter = n;
	while(true)
	{
		// Get the dual variables
		IloNumArray pi(env, n);
		rmp.getDuals(pi, partition_constraint);

		// for (size_t i = 0; i < n; i++)
		// {
		// 	cout << "Dual variable of constraint " << i << " = " << pi[i] << endl;
		// }
		
		// Build and solve the pricing problem
		double pricing_obj_val = 0;

		IloModel pricing_model(env);
		IloBoolVarArray x(env, n);
		IloExpr pricing_obj(env);
		IloExpr pricing_weight(env);
		IloRangeArray capacity_constraint(env);
		IloCplex pricing_problem;
		if (use_mip) {
			pricing_obj += 1;
			for (int i = 0; i < n; i++)
			{
				char var_name[50];
				sprintf(var_name, "x%d", i);
	
				x[i].setName(var_name);
				pricing_obj += -pi[i] * x[i];
				pricing_weight += weight[i] * x[i];
			}
			capacity_constraint.add(pricing_weight <= capacity);
	
			pricing_model.add(capacity_constraint);
	
			IloObjective pricing_objective = IloMinimize(env, pricing_obj);
			pricing_model.add(pricing_objective);
	
			pricing_problem = IloCplex(pricing_model);
			pricing_problem.setOut(env.getNullStream()); // disables CPLEX log
			pricing_problem.solve();
	
			pricing_obj_val = pricing_problem.getObjValue();
		}
		else {
			for (int i = 0; i < n; i++)
			{
				items[i].w = weight[i];
				items[i].p = max(pi[i],0.0)*M;
				items[i].x = 0;
				items[i].index = i;
			}

			pricing_obj_val = 1-(combo(&items[0], &items[n-1], capacity, 0, 0, true, false)/M);
		}
		
		cout << "Reduced cost is equal to " << pricing_obj_val << endl;
		if (pricing_obj_val < -1e-5)
		{

			// cout << "Reduced cost is equal to " << pricing_obj_val << ", which is less than 0..." << endl;

			IloNumArray entering_col(env, n);

			if (use_mip)
				pricing_problem.getValues(x, entering_col);
			else {
				for (size_t i = 0; i < n; i++)
				{
					if (items[i].x) {
						entering_col[items[i].index] = 1;
					}
					else
						entering_col[items[i].index] = 0;
				}
			}

			columns.push_back(vector<bool>(n, 0));
			cout << endl << "Entering column:" << endl;
			for (size_t i = 0; i < n; i++)
			{
				cout << (entering_col[i] < 0.5 ? 0 : 1) << " ";
				columns.back()[i] = (entering_col[i] < 0.5 ? 0 : 1);
			}
			cout << endl;

			// Add the column to the master problem
			// (the cost of the new variable is always 1)
			char var_name[50];
			sprintf(var_name, "y%d", lambda_counter++);
			IloNumVar new_lambda(master_objective(1) + partition_constraint(entering_col), 0, IloInfinity);
			new_lambda.setName(var_name);

			lambda.add(new_lambda);

			cout << "Solving the RMP again..." << endl;

			rmp.solve();

			cout << "New lower bound: " << rmp.getObjValue() << endl;

			cout << "New solution: ";
			for (size_t j = 0; j < lambda.getSize(); j++)
			{
				cout << rmp.getValue(lambda[j]) << " ";
			}
			cout << endl;

			pricing_problem.end();
			pricing_model.end();
			capacity_constraint.end();
			pricing_weight.end();
			pricing_obj.end();
			x.end();
		}
		else
		{
			pricing_problem.end();
			pricing_model.end();
			capacity_constraint.end();
			pricing_weight.end();
			pricing_obj.end();
			x.end();
			// cout << "No column with negative reduced costs found. The current basis is optimal" << endl;
			// cout << "Final master problem: " << endl;
			// system("cat model.lp");
			break;
		}
	}

	cout << endl;

	IloNumArray coeffs = IloNumArray(env, lambda.getSize());
	for (int i = 0; i < lambda.getSize(); i++)
	{
		coeffs[i] = 1;
	}
	master_objective.setLinearCoefs(lambda, coeffs);

	rmp.solve();

	vector<vector<double>> pairs(n, vector<double>(n,0));
	size_t fraci, fracj;
	bool all_int = AllInt(rmp, lambda, columns, n, pairs, fraci, fracj);
	cout << rmp.getObjValue() << " " << fraci << "," << fracj << endl;
	size_t used_lambdas = 0;

	list<Node> tree;

	if (!all_int) {
		Node n1;
		n1.used_columns = unordered_set<size_t>();
		n1.pairs = vector<pair<pair<size_t, size_t>, bool>>();
		n1.pairs.push_back({{fraci, fracj}, false});
		n1.LB = rmp.getObjValue();
		
		Node n2;
		n2.used_columns = unordered_set<size_t>();
		n2.pairs = vector<pair<pair<size_t, size_t>, bool>>();
		n2.pairs.push_back({{fraci, fracj}, true});
		n2.LB = rmp.getObjValue();
		
		for (size_t col = 0; col < columns.size(); col++) {
			n1.used_columns.insert(col);
			n2.used_columns.insert(col);
		}

		tree.push_back(n1);
		tree.push_back(n2);
	}

	while (!tree.empty()) {
		cout << "New node: ";
		Node nd = tree.back();
		tree.pop_back();
		for (auto pair = nd.pairs.cbegin(); pair != nd.pairs.cend(); pair++) {
			cout << pair->second;
		}
		cout << "\t";

		if (ceil(nd.LB) >= UB)
			continue;

		// for (auto it = nd.used_columns.cbegin(); it != nd.used_columns.cend(); it++)
		// 	cout << *it << " ";
		// cout << "\n";

		used_lambdas = 0;

		bool left_right = nd.pairs.back().second;
		fraci = nd.pairs.back().first.first; fracj = nd.pairs.back().first.second;
		for (size_t l = 0; l < lambda.getSize(); l++) {
			lambda[l].setUB(IloInfinity);


			if (left_right && (columns[l][fraci] != columns[l][fracj]))
				nd.used_columns.erase(l);
			else if (!left_right && columns[l][fraci] && columns[l][fracj])
				nd.used_columns.erase(l);

			if (nd.used_columns.find(l) == nd.used_columns.end())
				lambda[l].setUB(0);
			else
				used_lambdas++;
			
			// cout << lambda[l].getUB() << " ";
		}
		cout << used_lambdas << "\t";
		
		while(rmp.solve())
		{
			// Get the dual variables
			IloNumArray pi(env, n);
			rmp.getDuals(pi, partition_constraint);
			
			// for (size_t i = 0; i < n && nd.pairs.size() == 1; i++)
			// {
			// 	cout << "Dual variable of constraint " << i << " = " << pi[i] << "\n";
			// }
			
			// Build and solve the pricing problem
			double pricing_obj_val = 0;

			IloModel pricing_model(env);
			IloBoolVarArray x(env, n);
			IloExpr pricing_obj(env);
			IloExpr pricing_weight(env);
			IloRangeArray capacity_constraint(env);
			IloCplex pricing_problem;

			pricing_obj += 1;
			for (int i = 0; i < n; i++)
			{
				char var_name[50];
				sprintf(var_name, "x%d", i);
	
				x[i].setName(var_name);
				pricing_obj += -pi[i] * x[i];
				pricing_weight += weight[i] * x[i];
			}
			capacity_constraint.add(pricing_weight <= capacity);
	
			pricing_model.add(capacity_constraint);

			for (auto pair = nd.pairs.cbegin(); pair != nd.pairs.cend(); pair++) {
				if(pair->second) {
					pricing_model.add(x[pair->first.first] == x[pair->first.second]);
				}
				else {
					pricing_model.add(x[pair->first.first] + x[pair->first.second] <= 1);
				}
			}
	
			IloObjective pricing_objective = IloMinimize(env, pricing_obj);
			pricing_model.add(pricing_objective);
	
			pricing_problem = IloCplex(pricing_model);
			pricing_problem.setOut(env.getNullStream()); // disables CPLEX log
			pricing_problem.solve();
	
			pricing_obj_val = pricing_problem.getObjValue();
			
			// cout << "Reduced cost is equal to " << pricing_obj_val << "\n";
			if (pricing_obj_val < -1e-5)
			{
				IloNumArray entering_col(env, n);

				pricing_problem.getValues(x, entering_col);

				columns.push_back(vector<bool>(n, 0));

				nd.used_columns.insert(columns.size() - 1);
				// cout << "Entering column:\n";
				for (size_t i = 0; i < n; i++)
				{
					// cout << (entering_col[i] < 0.5 ? 0 : 1) << " ";
					columns.back()[i] = (entering_col[i] < 0.5 ? 0 : 1);
				}
				// cout << "\n";

				// Add the column to the master problem
				// (the cost of the new variable is always 1)
				char var_name[50];
				sprintf(var_name, "y%d", lambda_counter++);
				IloNumVar new_lambda(master_objective(1) + partition_constraint(entering_col), 0, IloInfinity);
				new_lambda.setName(var_name);

				lambda.add(new_lambda);

				pricing_problem.end();
				pricing_model.end();
				capacity_constraint.end();
				pricing_weight.end();
				pricing_obj.end();
				x.end();
			}
			else
			{
				pricing_problem.end();
				pricing_model.end();
				capacity_constraint.end();
				pricing_weight.end();
				pricing_obj.end();
				x.end();
				break;
			}
		}

		if (rmp.getStatus() == IloAlgorithm::Status::Infeasible) {
			cout << "\n";
			continue;
		}

		all_int = AllInt(rmp, lambda, columns, n, pairs, fraci, fracj);
		cout << fraci << "," << fracj << "\n";

		if (!all_int) {
			Node n1;
			n1.used_columns = nd.used_columns;
			n1.pairs = nd.pairs;
			n1.pairs.push_back({{fraci, fracj}, false});
			n1.LB = rmp.getObjValue();
			tree.push_back(n1);

			Node n2;
			n2.used_columns = nd.used_columns;
			n2.pairs = nd.pairs;
			n2.pairs.push_back({{fraci, fracj}, true});
			n2.LB = rmp.getObjValue();
			tree.push_back(n2);
		}
		else {
			UB = min(UB, (long long int)round(rmp.getObjValue()));
			cout << "Found int solution: " << rmp.getObjValue() << "\n";
		}
		cout << "\n";
	}

	for (size_t k = 0; k < lambda.getSize(); k++)
	{
		for (size_t j = 0; j < n; j++) {
			cout << columns[k][j] << " ";
		}
		cout << rmp.getValue(lambda[k]) << "\n";
	}
	cout << "Final solution: " << rmp.getObjValue() << "\n";
	
	env.end();

	return 0;
}

/*
- colocar combo.c
	- está dando diferente do MIP, tenho que olhar direito depois
- branching:
	- calcular valores de x baseados em lambdas
	- selecionar o par de x mais fracionário
	- à esquerda, esses x estão separados (lambdas com eles juntos = 0), add restrição no subproblema
		- soma dos x tem que dar <= 1
	- à direita, esses x estão juntos (lambdas com só 1 deles = 0), add restrição no subproblema
		- x = outro x
*/
