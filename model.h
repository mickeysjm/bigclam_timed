#ifndef MODEL_H
#define MODEL_H
#include <vector>
#include "graph.h"
#include "utils.h"
using namespace std;


struct Label {
	int id, date;
	string name;
	bool operator<(const Label &label2) {
		return date < label2.date;
	}
};

class Model {
public:
	Model() = default;
	~Model() = default;
	void Init(Utils & utils);
	void LoadGraph();
	void LoadLabel();
	int DetectCommunityCount();
	int GetOptimalCommunityCount();
	int CountDifferent(vector<int> v1, vector<int> v2);
	void Resize();
	void InitNeighbourCommunity(int community_count);
	void CalculateEtaF(vector<vector<double>> &eta);
	void AccumulateSumF();
	void AccumulateSumEtaF();
	void BackAccumulateSumF();
	double Norm2(vector<double> &v);
	double Norm2(vector<vector<double>> &m);
	void FGradientForRow(int u, vector<double> &gradient);
	double FLikelihoodForRow(int u, vector<double> &Fu);
	bool FGoStepSizeByLineSearch(int u, vector<double> &gradient, const int MAX_ITER = 10);
	void UpdateF(vector<int> &order);
	void EtaGradient(vector<vector<double>> &gradient);
	double EtaLikelihood(vector<vector<double>> &eta);
	bool EtaGoStepSizeByLineSearch(vector<vector<double>> &gradient, const int MAX_ITER = 10);
	void UpdateEta();
	void MLEGradAscent(const int MAX_ITER = 100, const int RECORD_PERIOD = 10);
	string GenerateModelName(int iteration);
	void SaveF(string &model_name);
	void SaveEta(string &model_name);
	void SaveData(string model_name);

private:
	//parameters
	string output_file_perfix_;	// Output Graph data prefix
	string input_graph_filename_;	// Input edgelist file name
	string input_label_filename_;	// Input file name for node names (Node ID, Node label)
	int optimal_community_count_;	// The number of communities to detect (-1: detect automatically)
	int min_community_count_;	// Minimum number of communities to try
	int max_community_count_;	// Maximum number of communities to try
	int divide_community_count_;	// How many trials for the number of communities
	int thread_number_;	// Number of threads for parallelization
	double step_alpha_;	// Alpha for backtracking line search
	double step_beta_;	// Beta for backtracking line search
	double MIN_F_; // max value for F
	double MAX_F_; // min value for F
	double MIN_ETA_; // max value for Eta
	double MAX_ETA_; // min value for Eta
	double MIN_P_; // min probability for P(u, v)
	double MAX_P_; // max probability for P(u, v)

	int community_count_; // Number of Community
	Graph graph_;	// Input graph data
	vector<vector<double>> F_; // Strength of affiliation edge between community and node, Nodes * Coms
	vector<vector<double>> newF_; // new F_, Nodes * Coms
	vector<vector<double>> Eta_; // Coms * Coms
	vector<vector<double>> newEta_; // Coms * Coms
	vector<vector<double>> EtaF_; // eta times Fv, Nodes * Coms
	vector<Label> label_; // message about each node, Nodes
	map<int, vector<double>> sumF_; // accumulated sum of F, Dates * Coms
	map<int, vector<double>> sumEtaF_; // accumulated sum of EtaF, Dates * Coms
	map<int, vector<double>> backsumF_; // accumulated sum from tail to head of F, Dates * Coms
	map<int, int> id2order_; // convert paper id to order.
};
#endif
