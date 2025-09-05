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
                cerr << k << "/" << n << " break " << cut_val << "\n";
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

            cut_val += 2 -2 * maxback_val[i]; // TODO If this is 0, algorithm can be stopped because a subtour was found.
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
    return vector<vector<int>>();
}
