#include "Kruskal.h"

Kruskal::Kruskal(const vvi& dist, int zero_setup) {
	this->zero_setup = zero_setup;
	for (size_t i = (zero_setup == USE_ZERO ? 0 : 1); i < dist.size(); ++i) {
		for (size_t j = i+1; j < dist[i].size(); ++j) {
			graph.push( make_pair(-dist[i][j], make_pair(i, j)) );
		}
	}
	zero_edges = vector<double>(dist.size());
	for(size_t i = 0; i < dist.size(); ++i) {
		zero_edges[i] = dist[0][i];
	}
}

void Kruskal::initDisjoint(size_t n){
	pset.resize(n);
	for (size_t i = 0; i < n; ++i) {
		pset[i] = i;
	}
}

size_t Kruskal::findSet(size_t i){
	return (pset[i] == i) ? i : (pset[i] = findSet(pset[i]));
}

void Kruskal::unionSet(size_t i, size_t j) {
	pset[findSet(i)] = findSet(j);
}

bool Kruskal::isSameSet(size_t i, size_t j){
	return (findSet(i) == findSet(j))? true:false;
}

vector<pair<size_t,size_t>> Kruskal::getEdges(){
	return edges;
}

double Kruskal::MST(size_t nodes){
	initDisjoint(nodes);

	double cost = 0;

	while(!graph.empty()){
		pair<double, ii> p = graph.top();
		graph.pop();

		if(!isSameSet(p.second.first, p.second.second)){
			edges.push_back(make_pair(p.second.first, p.second.second));
			cost += (-p.first);
			unionSet(p.second.first, p.second.second);
		}
	}
	if (zero_setup == TWO_CLOSEST) {
		size_t closest = 1;
		size_t second_closest = 2;
		for (size_t i=3; i<zero_edges.size(); i++) {
			if (zero_edges[i] < zero_edges[second_closest]) {
				if (zero_edges[i] < zero_edges[closest]) {
					second_closest = closest;
					closest = i;
				}
				else
					second_closest = i;
			}
		}
		edges.push_back(make_pair(0, closest));
		cost += (zero_edges[closest]);
		unionSet(0, closest);
		edges.push_back(make_pair(0, second_closest));
		cost += (zero_edges[second_closest]);
		unionSet(0, second_closest);
	}

	return cost;
}