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

vector <vector<int> > MinCut(double** x, int n) {
    // cerr << "mincut\n";
    vector<vector<int>> Ss = vector<vector<int>>(2, vector<int>());
    double mincut_val = MAXFLOAT;
    vector<int> disjointed_set(n);
    double** xx = new double*[n];
	for (int i = 0; i < n; i++) {
        disjointed_set[i] = i;
		xx[i] = new double[n];
        for (int j = i; j < n; j++) {
            xx[i][j] = x[i][j];
        }
	}

    vector<double> degree(n,0);

    int nn = n;

    for (int k = 1; k < n; k++) {
        //Below: cut_val, s, t = MaxBackStep();
        vector<double> maxback_val = vector<double>(nn);
        vector<int> insertion_order(nn,0);
        vector<bool> Sk(nn, false);
        Sk[0] = true;

        for (int i = 0; i < nn; i++) {
            maxback_val[i] = xx[0][i];
            degree[i] = 0;
        }
        
        for (int i = 0; i < nn; i++) {
            for (int j = i+1; j < nn; j++) {
                degree[i] += xx[i][j];
                degree[j] += xx[i][j];
            }
        }

        double cut_val = 0;
        for (int i = 1; i < nn; i++) {
            cut_val += xx[0][i];
        }

        for (int kk = 1; kk < nn; kk++) {
            int i = -1;
            bool first = true;
            for (int j = 1; j < nn; j++) { // we can ignore 0 because kk starts from 0, thus Sk[0] will always be true.
                if (Sk[j])
                    continue;
                
                if (first || maxback_val[j] > maxback_val[i]) {
                    i = j;
                    first = false;
                }
            }
            Sk[i] = true;
            insertion_order[kk] = i;

            cut_val += degree[i] -2 * maxback_val[i];
            for (int j = 0; j < nn; j++) {
                if (Sk[j])
                    continue;
                maxback_val[j] += (i < j) ? (xx[i][j]) : (xx[j][i]);
            }
        }
        
        int s = (*(insertion_order.rbegin()+1));
        int t = (*insertion_order.rbegin());

        if (cut_val + EPSILON < mincut_val) {
            mincut_val = cut_val;
            Ss[0] = vector<int>();
            for (int j = 0; j < n; j++) {
                if (disjointed_set[j] == t)
                    Ss[0].push_back(j);
            }
        }

        if (t < s)
            swap(s, t);
        
        for (int j = 0; j < n; j++)
            if (disjointed_set[j] == t)
                disjointed_set[j] = s;

        nn--;
        double** xxx = new double*[nn];
        for (int i = 0; i < nn; i++) {
            xxx[i] = new double[nn];
            for (int j = i; j < nn; j++) {
                if (i != s && i != t && j != s && j != t)
                    xxx[i][j] = xx[i][j];
                else {
                    xxx[i][j] = (s<j ? xx[s][j] : xx[j][s]) + (t<j ? xx[t][j] : xx[j][t]);
                }
            }
        }

        for (int i = 0; i < nn+1; i++) {
            delete[] xx[i];
        }
        delete[] xx;

        xx = xxx;
        xxx = NULL;

        if (abs(mincut_val) < EPSILON) {
            Ss = vector<vector<int>>();
            break;
        }
    }

    // if (Ss[0].size() == 0 || Ss[0].size() == n)
    //     Ss[0] = {0};

    for (int i = 0; i < nn; i++) {
		delete[] xx[i];
	}
	delete[] xx;

    if (Ss.size() == 0)
        return Ss;
    
    vector<bool> inSs(n, false);
    for (size_t k = 0; k < Ss[0].size(); k++) {
        inSs[Ss[0][k]] = true;
    }

    for (int i = 0; i < n; i++)
        if (!inSs[i])
            Ss[1].push_back(i);

    return Ss;
}
