#ifndef GRAPH_H
#define GRAPH_H
#include <vector>
#include <cassert>
using namespace std;


class Graph {
public:
	Graph() = default;
	Graph(int node_count) : node_count_(node_count), out_edges_(node_count), in_edges_(node_count) {}
	~Graph() = default;
	void AddEdge(int from, int to);
	void SetNodeCount(int node_count);
	int GetNodeCount(){
		return node_count_;
	}
	vector<int> &GetInNeighbour(int index){
		assert((0 <= index && index < in_edges_.size()) && "index out of range");
		return in_edges_[index];
	}
	vector<int> &GetOutNeighbour(int index){
		assert((0 <= index && index < out_edges_.size()) && "index out of range");
		return out_edges_[index];
	}
private:
	vector<vector<int>> out_edges_, in_edges_;
	int node_count_;
};
#endif