#include "Kruskal.h"

Kruskal::Kruskal(const vvi& _dist, int zero_setup, vector <size_t>& _pset, vector<ii>& _non_zero_edges, vector<pair<size_t,size_t>>* _edges) {
	dist = &_dist;
	pset = &_pset;
	non_zero_edges = &_non_zero_edges;
	edges = _edges;
	// (*non_zero_edges).reserve(_dist.size() * _dist.size() / 2);
	// this->zero_setup = zero_setup;
	// for (size_t i = (zero_setup == USE_ZERO ? 0 : 1); i < _dist.size(); ++i) {
	// 	for (size_t j = i+1; j < _dist[i].size(); ++j) {
	// 		if (_dist[i][j] < 9e7)
	// 			non_zero_edges.push_back( make_pair(i, j) );
	// 	}
	// }
	sort((*non_zero_edges).begin(), (*non_zero_edges).end(), [&](const ii& a, const ii& b){ return _dist[a.first][a.second] < _dist[b.first][b.second]; });
}

Kruskal::~Kruskal(){
	non_zero_edges = NULL;
	pset = NULL;
	dist = NULL;
	edges = NULL;
}

void Kruskal::initDisjoint(size_t n){
	for (size_t i = 0; i < n; ++i) {
		(*pset)[i] = i;
	}
}

size_t Kruskal::findSet(size_t i){
	return ((*pset)[i] == i) ? i : ((*pset)[i] = findSet((*pset)[i]));
}

void Kruskal::unionSet(size_t i, size_t j) {
	(*pset)[findSet(i)] = findSet(j);
}

bool Kruskal::isSameSet(size_t i, size_t j){
	return (findSet(i) == findSet(j))? true:false;
}

double Kruskal::MST(size_t nodes){
	initDisjoint(nodes);

	double cost = 0;
	size_t graphsize = 0;

	for (auto it = (*non_zero_edges).cbegin(); it != (*non_zero_edges).cend(); it++) {
		if(!isSameSet(it->first, it->second)) {
			(*edges)[graphsize] = make_pair(it->first, it->second);
			graphsize++;
			cost += (*dist)[it->first][it->second];
			unionSet(it->first, it->second);
		}
	}

	size_t closest = 1;
	size_t second_closest = 2;
	if ((*dist)[0][2] < (*dist)[0][1])
		swap(closest, second_closest);
	for (size_t i=3; i<dist[0].size(); i++) {
		if ((*dist)[0][i] < (*dist)[0][second_closest]) {
			if ((*dist)[0][i] < (*dist)[0][closest]) {
				second_closest = closest;
				closest = i;
			}
			else
				second_closest = i;
		}
	}
	(*edges)[graphsize] = make_pair(0, closest);
	graphsize++;
	cost += ((*dist)[0][closest]);
	// unionSet(0, closest);
	(*edges)[graphsize] = make_pair(0, second_closest);
	graphsize++;
	cost += ((*dist)[0][second_closest]);
	// unionSet(0, second_closest);
	// if (graphsize != (*edges).size())
	// 	cerr << graphsize << " diff\n";

	return cost;
}