#include <ilcplex/ilocplex.h>
#include "combo.c"
#include <vector>
#include <cstdlib>

using namespace std;

int main() 
{
	const double M = 1e6;
	vector<int> weight = {2, 1, 3, 3, 5};
	int capacity = 7;
	int n = weight.size();
	bool combomip = false;

	item* items = new item[weight.size()];
	for (int i = 0; i < n; i++)
	{
		items[i].w = weight[i];
		items[i].profitofitem = M;
		items[i].x = 0;
		items[i].index = i;
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

		for (size_t i = 0; i < n; i++)
		{
			cout << "Dual variable of constraint " << i << " = " << pi[i] << endl;
		}
		
		// Build and solve the pricing problem
		double pricing_obj_val = 0;

		IloModel pricing_model(env);
		IloBoolVarArray x(env, n);
		IloExpr pricing_obj(env);
		IloExpr pricing_weight(env);
		IloRangeArray capacity_constraint(env);
		IloCplex pricing_problem;
		if (combomip) {
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
				items[i].profitofitem = max(pi[i],0.0);
				items[i].x = 0;
				items[i].index = i;
			}

			pricing_obj_val = 1-combo(&items[0], &items[n-1], capacity, 0, 0, true, false);
		}
		
		cout << "Reduced cost is equal to " << pricing_obj_val << endl;
		if (pricing_obj_val < -1e-5)
		{

			// cout << "Reduced cost is equal to " << pricing_obj_val << ", which is less than 0..." << endl;

			IloNumArray entering_col(env, n);

			if (combomip)
				pricing_problem.getValues(x, entering_col);

			cout << endl << "Entering column:" << endl;
			if (!combomip) {
				for (size_t i = 0; i < n; i++)
				{
					if (items[i].x) {
						entering_col[items[i].index] = 1;
					}
					else
						entering_col[items[i].index] = 0;
				}
			}
			for (size_t i = 0; i < n; i++)
			{
				cout << (entering_col[i] < 0.5 ? 0 : 1) << endl;
			}
			cout << endl;

			// Add the column to the master problem
			// (the cost of the new variable is always 1)
			char var_name[50];
			sprintf(var_name, "y%d", lambda_counter++);
			IloNumVar new_lambda(master_objective(1) + partition_constraint(entering_col), 0, IloInfinity);
			// IloNumVar new_lambda2(master_objective(1) + partition_constraint(entering_col), 0, IloInfinity);
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

			lambda_counter++;

		}
		else
		{
			// cout << "No column with negative reduced costs found. The current basis is optimal" << endl;
			cout << "Final master problem: " << endl;
			system("cat model.lp");
			break;
		}
	}

	cout << endl;
	cout << "Forcing items 1 and 2 to be separated in the master (for branch-and-price only): " << endl;
	// 0 1 2 3 4 5 6 7 8 9 10 11
	//                           
	// 1 0 0 0 0 1 1 1 0 1  0  0
	// 0 1 0 0 0 1 1 0 0 0  1  1
	// 0 0 1 0 0 1 0 1 1 0  0  1
	// 0 0 0 1 0 0 1 0 1 0  0  1
	// 0 0 0 0 1 0 0 0 0 1  1  0
	//           v v v f v  f  f


	// itens 1 and 2 are together only on columns 5 and 11
	lambda[11].setUB(0.0);
	lambda[5].setUB(0.0);

	// to allow them again:
	// lambda[5].setUB(IloInfinity);
	// lambda[11].setUB(IloInfinity);

	rmp.solve();

	for (size_t j = 0; j < lambda.getSize(); j++)
	{
		cout << rmp.getValue(lambda[j]) << " ";
	}
	cout << endl;
	env.end();

	return 0;
}
/*
- colocar combo.c
- branching:
	- calcular valores de x baseados em lambdas
	- selecionar os 2 x mais fracionários
	- à esquerda, esses x estão separados (lambdas com eles juntos = 0), add restrição no subproblema
	- à direita, esses x estão juntos (lambdas com só 1 deles = 0), add restrição no subproblema
*/