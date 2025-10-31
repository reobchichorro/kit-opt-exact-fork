#include "dataFunction.h"
#include "auxFunctions.h"
#include "separation.h"
#include "MyLazyCallback.h"
#include "MyCutCallback.h"
#include "MyBranchCallback.h"
#include <stdio.h>
#include <iostream>
#include <unordered_map>
#include <ilcplex/ilocplex.h>

void setUBs(unordered_map<string, double>& UBs) {
	UBs["a280"] = 2579;
	UBs["ali535"] = 202339;
	UBs["att48"] = 10628;
	UBs["att532"] = 27686;
	UBs["bayg29"] = 1610;
	UBs["bays29"] = 2020;
	UBs["berlin52"] = 7542;
	UBs["bier127"] = 118282;
	UBs["brazil58"] = 25395;
	UBs["brd14051"] = 469385;
	UBs["brg180"] = 1950;
	UBs["burma14"] = 3323;
	UBs["ch130"] = 6110;
	UBs["ch150"] = 6528;
	UBs["d198"] = 15780;
	UBs["d493"] = 35002;
	UBs["d657"] = 48912;
	UBs["d1291"] = 50801;
	UBs["d1655"] = 62128;
	UBs["d2103"] = 80450;
	UBs["d15112"] = 1573084;
	UBs["d18512"] = 645238;
	UBs["dantzig42"] = 699;
	UBs["dsj1000euc_2d"] = 18659688;
	UBs["dsj1000ceil_2d"] = 18660188;
	UBs["eil51"] = 426;
	UBs["eil76"] = 538;
	UBs["eil101"] = 629;
	UBs["fl417"] = 11861;
	UBs["fl1400"] = 20127;
	UBs["fl1577"] = 22249;
	UBs["fl3795"] = 28772;
	UBs["fnl4461"] = 182566;
	UBs["fri26"] = 937;
	UBs["gil262"] = 2378;
	UBs["gr17"] = 2085;
	UBs["gr21"] = 2707;
	UBs["gr24"] = 1272;
	UBs["gr48"] = 5046;
	UBs["gr96"] = 55209;
	UBs["gr120"] = 6942;
	UBs["gr137"] = 69853;
	UBs["gr202"] = 40160;
	UBs["gr229"] = 134602;
	UBs["gr431"] = 171414;
	UBs["gr666"] = 294358;
	UBs["hk48"] = 11461;
	UBs["kroA100"] = 21282;
	UBs["kroB100"] = 22141;
	UBs["kroC100"] = 20749;
	UBs["kroD100"] = 21294;
	UBs["kroE100"] = 22068;
	UBs["kroA150"] = 26524;
	UBs["kroB150"] = 26130;
	UBs["kroA200"] = 29368;
	UBs["kroB200"] = 29437;
	UBs["lin105"] = 14379;
	UBs["lin318"] = 42029;
	UBs["linhp318"] = 41345;
	UBs["nrw1379"] = 56638;
	UBs["p654"] = 34643;
	UBs["pa561"] = 2763;
	UBs["pcb442"] = 50778;
	UBs["pcb1173"] = 56892;
	UBs["pcb3038"] = 137694;
	UBs["pla7397"] = 23260728;
	UBs["pla33810"] = 66048945;
	UBs["pla85900"] = 142382641;
	UBs["pr76"] = 108159;
	UBs["pr107"] = 44303;
	UBs["pr124"] = 59030;
	UBs["pr136"] = 96772;
	UBs["pr144"] = 58537;
	UBs["pr152"] = 73682;
	UBs["pr226"] = 80369;
	UBs["pr264"] = 49135;
	UBs["pr299"] = 48191;
	UBs["pr439"] = 107217;
	UBs["pr1002"] = 259045;
	UBs["pr2392"] = 378032;
	UBs["rat99"] = 1211;
	UBs["rat195"] = 2323;
	UBs["rat575"] = 6773;
	UBs["rat783"] = 8806;
	UBs["rd100"] = 7910;
	UBs["rd400"] = 15281;
	UBs["rl1304"] = 252948;
	UBs["rl1323"] = 270199;
	UBs["rl1889"] = 316536;
	UBs["rl5915"] = 565530;
	UBs["rl5934"] = 556045;
	UBs["rl11849"] = 923288;
	UBs["si175"] = 21407;
	UBs["si535"] = 48450;
	UBs["si1032"] = 92650;
	UBs["st70"] = 675;
	UBs["swiss42"] = 1273;
	UBs["ts225"] = 126643;
	UBs["tsp225"] = 3916;
	UBs["u159"] = 42080;
	UBs["u574"] = 36905;
	UBs["u724"] = 41910;
	UBs["u1060"] = 224094;
	UBs["u1432"] = 152970;
	UBs["u1817"] = 57201;
	UBs["u2152"] = 64253;
	UBs["u2319"] = 234256;
	UBs["ulysses16"] = 6859;
	UBs["ulysses22"] = 7013;
	UBs["usa13509"] = 19982859;
	UBs["vm1084"] = 239297;
	UBs["vm1748"] = 336556;
	for (auto it = UBs.begin(); it != UBs.end(); it++)
		UBs[it->first] = it->second + 1;
}

void STSP_Solve(Data *data, string instanceName, double ub)
{   
    IloEnv env;
    IloModel model(env);

    env.setName("Branch and Cut");
    model.setName("Symmetrical Traveling Salesman Problem");

    int dimension = data->getDimension();

    /********** Creating variable x for each edge **********/
    IloArray <IloBoolVarArray> x(env, dimension);

    for (int i = 0; i < dimension; i++) {
		IloBoolVarArray array(env, dimension);
		x[i] = array;
	}
    /*******************************************************/

    /*********** Adding x variables to the model ***********/
    char var[100];
    for (int i = 0; i < dimension; i++){
        for (int j = i + 1; j < dimension; j++){
            sprintf(var, "X(%d,%d)", i, j);
			x[i][j].setName(var);
			model.add(x[i][j]);
        }
    }
    /******************************************************/
    
    /**************** Objective Function ******************/
    IloExpr obj(env);
    for (int i = 0; i < dimension; i++) {	
		for (int j = i + 1; j < dimension; j++) {
			obj += data->getDistance(i, j)*x[i][j];
		}
	}
    model.add(IloMinimize(env, obj));
    /******************************************************/
    
    /******************** Constraints *********************/
    IloRange r;
    char name[100];

    for (int i = 0; i < dimension; i++){
        IloExpr sumX(env);
        for (int j = 0; j < dimension; j++){
            if (j < i) {
				sumX += x[j][i];
			}
            if (i < j){
                sumX += x[i][j];
            }
        }
        r = (sumX == 2);
        sprintf(name, "c_%d", i);
		r.setName(name);
		model.add(r);
    }
    /******************************************************/

    /****************** Solve the model *******************/
    IloCplex STSP(model);
    STSP.setParam(IloCplex::TiLim, 2*60*60);
    STSP.setParam(IloCplex::Threads, 1);
    STSP.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 1e-08);
    STSP.setParam(IloCplex::CutUp, ub);
	// STSP.setParam(IloCplex::Param::MIP::Strategy::VariableSelect, CPX_VARSEL_STRONG);
    //STSP.exportModel("stsp.lp");

    double timeBefore, timeAfter;

    const IloArray<IloBoolVarArray>& x_ref = x;

    /********** Creating Branch Callback Object ***********/
    MyBranchCallback* branchCbk = new (env) MyBranchCallback(env);
    STSP.use(branchCbk);
    /******************************************************/

    /************ Creating Cut Callback Object ************/
    MyCutCallback* cutCbk = new (env) MyCutCallback(env, x_ref, dimension); 
	STSP.use(cutCbk);
    /******************************************************/

    /************ Creating Lazy Callback Object ***********/
    MyLazyCallback* lazyCbk = new (env) MyLazyCallback(env, x_ref, dimension);
    STSP.use(lazyCbk);
    /******************************************************/

    try{ 
	    timeBefore = STSP.getTime();
	    STSP.solve();
	    timeAfter = STSP.getTime();
    }
    catch(IloException& e){
        std::cout << e;
    }

    printResults(STSP, instanceName, timeAfter-timeBefore);
	printSolution(STSP, x, dimension);
	//printResultsToFile(STSP, instanceName, timeAfter-timeBefore);
    /******************************************************/

    /**************** Cleaning the memory *****************/
    delete branchCbk;
    delete cutCbk;
    delete lazyCbk;
    env.end();
    /******************************************************/
}

int main(int argc, char** argv)
{
    double ub = numeric_limits<double>::max();

    if(argc < 2){
        printf("Correct command: ./bc data/instance\n");
        return 0;
    }
    else if(argc > 2){
        ub = stod(argv[2]);
    }

    unordered_map<string, double> UBs;
	setUBs(UBs);

    Data * data = new Data(argc, argv[1]);
    data->readData();

    ub = UBs[data->getInstanceName()];

    string instanceName = getInstanceName(argv);
    STSP_Solve(data, instanceName, ub);

    /*************** Cleaning the memory ***************/
    delete data;
    /***************************************************/

    return 0;
}