#include "graph.h"

void Graph::AddEdge(int from, int to) {
	out_edges_[from].push_back(to);
	in_edges_[to].push_back(from);
}

void Graph::SetNodeCount(int node_count) {
	node_count_ = node_count;
	out_edges_.clear();
	out_edges_.resize(node_count_);
	in_edges_.clear();
	in_edges_.resize(node_count_);
}


// int main(int argc, char const *argv[])
// {
// 	Graph g(10);
// 	return 0;
// }