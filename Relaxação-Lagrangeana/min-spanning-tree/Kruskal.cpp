#include "Kruskal.h"

Kruskal::Kruskal(vvi dist){
	for(size_t i = 0; i < dist.size(); ++i){
		for(size_t j = i+1; j < dist[i].size(); ++j){
			graph.push( make_pair(-dist[i][j], make_pair(i, j)) );
		}	
	}
}

void Kruskal::initDisjoint(size_t n){
	pset.resize(n);
	for (size_t i = 0; i < n; ++i){
		pset[i] = i;
	}
}

size_t Kruskal::findSet(size_t i){
	return (pset[i] == i) ? i : (pset[i] = findSet(pset[i]));
}

void Kruskal::unionSet(size_t i, size_t j){
	pset[findSet(i)] = findSet(j);
}

bool Kruskal::isSameSet(size_t i, size_t j){
	return (findSet(i) == findSet(j))? true:false;
}

vii Kruskal::getEdges(){
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

	return cost;
}