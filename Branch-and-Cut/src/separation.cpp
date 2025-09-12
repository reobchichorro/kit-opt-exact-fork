#include "separation.h"

vector <vector<int> > MaxBack(double** x, int n) {
    // kit's MaxBack finds one subtour. The method here needs to return ALL subtours
    vector<vector<int>> Ss; // TODO find better name
    vector<bool> found = vector<bool>(n,false); // Which vertices are in one of Ss' vectors

    // To be reset at each iteration
    vector<double> maxback_val = vector<double>(n);

    // for (int i = 0; i < n; i++) {
    //     for (int j = i + 1; j < n; j++) {
    //         if (x[i][j] > 0.1)
    //             cerr << i << "-" << j << "," << x[i][j] << " ";
    //         else if (x[i][j] > EPSILON)
    //             cerr << "eps-" << i << "-" << j << "," << x[i][j] << " ";
    //     }
    // }
    // cerr << "\n";

    for (int s0 = 0; s0 < n; s0++) {
        if (found[s0])
            continue;
        found[s0] = true;
        vector<bool> Sk(n, false);
        Sk[s0] = true;
        
        for (int i = 0; i < s0; i++) {
            maxback_val[i] = x[i][s0];
        }
        for (int i = s0; i < n; i++) {
            maxback_val[i] = x[s0][i];
        }

        double cut_val = 0;
        for (int i = 0; i < s0; i++) {
            cut_val += x[i][s0];
        }
        for (int i = s0+1; i < n; i++) {
            cut_val += x[s0][i];
        }
        double mincut_val = cut_val;

        Ss.push_back({s0});

        for (int k = 1; k < n; k++) { // Always iterates n-1 times, which should be equal to |V|-|S_0|
            if (abs(cut_val) < EPSILON) {
                // cerr << k << "/" << n << " break " << cut_val << "\n";
                break;
            }

            int i = -1;
            bool first = true;
            for (int j = 0; j < n; j++) {
                if (Sk[j])
                    continue;
                
                if (first || maxback_val[j] > maxback_val[i]) {
                    i = j;
                    first = false;
                }
            }
            Sk[i] = true;

            cut_val += 2 -2 * maxback_val[i];
            for (int j = 0; j < n; j++) {
                if (Sk[j])
                    continue;
                maxback_val[j] += (i < j) ? (x[i][j]) : (x[j][i]);
            }

            if (cut_val + EPSILON < mincut_val) {
                mincut_val = cut_val;
                Ss.back() = vector<int>();
                for (int j = 0; j < n; j++)
                    if (Sk[j])
                        Ss.back().push_back(j);
            }
        }

        for (auto it = Ss.back().cbegin(); it != Ss.back().cend(); it++)
            found[*it] = true;
    }
    if (Ss.size() == 1 && Ss[0].size() == n)
        Ss = vector<vector<int>>();
        
    return Ss;
}

// void MinCutStep(double** x, int n, vector<int>& S, int s0) {
    
// }

vector <vector<int> > MinCut(double** x, int n) {
    // cerr << "mincut\n";
    vector<vector<int>> Ss = vector<vector<int>>();
    vector<int> disjointed_set(n);
    double** xx = new double*[n];
    list<int> remaining;

	for (int i = 0; i < n; i++) {
        disjointed_set[i] = i;
        remaining.push_back(i);
		xx[i] = new double[n];
        for (int j = i; j < n; j++) {
            xx[i][j] = x[i][j];
        }
	}

    // for (int i = 0; i < n; i++) {
    //     for (int j = i+1; j < n; j++) {
    //         if (x[i][j] > 0.1)
    //             cerr << i << "-(" << x[i][j] << ")-" << j << " ";
    //         else if (x[i][j] > EPSILON)
    //             cerr << "eps-" << i << "-(" << x[i][j] << ")-" << j << " ";
    //     }
    // }
    // cerr << "\n";

    vector<double> degree(n,0);

    int nn = n;
    int s0 = 0;

    for (int k = 1; k < n; k++) {
        //Below: cut_val, s, t = MaxBackStep();
        vector<double> maxback_val = vector<double>(n, -3); // Using -3 as a dummy value. Positions that have -3 should never be accessed.
        vector<int> insertion_order(nn,-1);
        vector<bool> Sk(n, false);
        Sk[s0] = true;
        insertion_order[0] = s0;

        for (auto it = remaining.cbegin(); it != remaining.cend(); it++) {
            if (s0 < *it) {
                maxback_val[*it] = xx[s0][*it];
            } else {
                maxback_val[*it] = xx[*it][s0];                
            }
            degree[*it] = 0;
        }

        for (auto it = remaining.cbegin(); it != remaining.cend(); it++) {
            for (auto itj = next(it,1); itj != remaining.cend(); itj++) {
                degree[*it] += xx[*it][*itj];
                degree[*itj] += xx[*it][*itj];
            }
        }

        double cut_val = 0;
        for (auto it = remaining.cbegin(); it != remaining.cend(); it++) {
            if (s0 < *it) {
                cut_val += xx[s0][*it];
            } else if (s0 > *it) {
                cut_val += xx[*it][s0];                
            }
        }

        for (int kk = 1; kk < nn; kk++) {
            int i = -1;
            bool first = true;
            
            for (auto it = remaining.cbegin(); it != remaining.cend(); it++) {
                if (Sk[*it])
                    continue;
                
                if (first || maxback_val[*it] > maxback_val[i]) {
                    i = *it;
                    first = false;
                }
            }
            Sk[i] = true;
            insertion_order[kk] = i;

            if (kk < nn-1)
                cut_val += degree[i] -2 * maxback_val[i];
            for (auto it = remaining.cbegin(); it != remaining.cend(); it++) {
                if (Sk[*it])
                    continue;
                maxback_val[*it] += (i < *it) ? (xx[i][*it]) : (xx[*it][i]);
            }
        }
        
        int s = (*(insertion_order.rbegin()+1));
        int t = (*insertion_order.rbegin());

        if (cut_val + EPSILON < 2) {
            size_t cut1 = Ss.size();
            Ss.push_back({});
            for (int j = 0; j < n; j++) {
                if (disjointed_set[j] == t)
                    Ss[cut1].push_back(j);
            }
            cerr << cut_val << ":\t";
            for (size_t j=0; j<Ss[cut1].size(); j++)
                cerr << Ss[cut1][j] << ",";
            cerr << "\n";
        }

        if (t < s)
            swap(s, t);
        
        for (int j = 0; j < n; j++)
            if (disjointed_set[j] == t)
                disjointed_set[j] = s;

        remaining.remove(t);

        nn--;
        double** xxx = new double*[n];
        for (int i = 0; i < n; i++) {
            xxx[i] = new double[n];
            for (int j = i; j < n; j++) {
                if (disjointed_set[i] != s) {
                    if (disjointed_set[j] != s)
                        xxx[i][j] = xx[i][j];
                    else
                        xxx[i][j] = (s<i ? xx[s][i] : xx[i][s]) + (t<i ? xx[t][i] : xx[i][t]);    
                }
                else {
                    xxx[i][j] = (s<j ? xx[s][j] : xx[j][s]) + (t<j ? xx[t][j] : xx[j][t]);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            delete[] xx[i];
        }
        delete[] xx;

        xx = xxx;
        xxx = NULL;
    }

    for (int i = 0; i < nn; i++) {
		delete[] xx[i];
	}
	delete[] xx;

    
    if (Ss.size() > 0) {
        vector<bool> inSs(n, false);
        for (size_t i = 0; i < Ss.size(); i++) {
            for (auto it = Ss[i].cbegin(); it != Ss[i].cend(); it++) {
                inSs[*it] = true;
            }
        }
        bool first = true;
        for (int i = 0; i < n; i++) {
            if (!inSs[i]) {
                if (first) {
                    first = false;
                    Ss.push_back(vector<int>());
                }
                Ss.back().push_back(i);
                cerr << "\t" << i << ",";
            }
            cerr << "\n";
        }
    }

    return Ss;
}
