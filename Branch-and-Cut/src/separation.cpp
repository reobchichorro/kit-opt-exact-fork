#include "separation.h"

vector <vector<int> > MaxBack(double** x, int n) {
    // kit's MaxBack finds one subtour. The method here needs to return ALL subtours
    vector<vector<int>> ans = vector<vector<int>>(1, vector<int>(1,0));
    vector<double> maxback_val = vector<double>(n,0);
    for (int i = 0; i < n; i++) {
        maxback_val[i] = x[0][i];
	}
    double cut_val = 0;
    for (int i = 1; i < n; i++) {
        cut_val += x[0][i];
	}
    double mincut_val = cut_val;

    for (int k = 1; k < n; k++) {
        int i = -1;
        double max_maxback_val = 0;
        for (int j = 1; j < n; j++) {
            bool skip = false;
            for (size_t jj = 0; jj < ans[k-1].size(); jj++) {
                if (j == ans[k-1][jj]) {
                    skip = true;
                    break;
                }
            }
            if (skip)
                continue;
            
            if (maxback_val[j] > max_maxback_val) {
                max_maxback_val = maxback_val[j];
                i = j;
            }
        }
        ans.push_back(vector<int>(ans[k-1].size()+1, 0));
        for (size_t ii = 0; i < ans[k-1].size(); ii++) {
            ans[k][ii] = ans[k-1][ii];
        }
        ans[k][ans[k-1].size()] = i;
        cut_val += 2 -2 * maxback_val[i]; // TODO substitute 2 for degree(i)

    }

    return ans;
}

vector <vector<int> > MinCut(double** x, int n) {
    return vector<vector<int>>();
}
