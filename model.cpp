#include "model.h"
#include "common.h"
#include <fstream>
#include <sstream>
#include <iterator>
#include <iostream>
#include <ctime>
#include <random>
#include <algorithm>
#include <iomanip>

void Model::Init(Utils & utils) {
	output_file_perfix_ = utils.ParseArgumentString("-o", "output\\saved_", "Output Graph data prefix");
	input_graph_filename_ = utils.ParseArgumentString("-i", "test_graph.in", "Input edgelist file name");
	input_label_filename_ = utils.ParseArgumentString("-l", "test_label.in", "Input file name for node names (Node ID, Node label) ");
	optimal_community_count_ = utils.ParseArgumentInt("-c", 100, "The number of communities to detect (-1 detect automatically)");
	min_community_count_ = utils.ParseArgumentInt("-mc", 5, "Minimum number of communities to try");
	max_community_count_ = utils.ParseArgumentInt("-xc", 100, "Maximum number of communities to try");
	divide_community_count_ = utils.ParseArgumentInt("-nc", 10, "How many trials for the number of communities");
	thread_number_ = utils.ParseArgumentInt("-nt", 4, "Number of threads for parallelization");
	step_alpha_ = utils.ParseArgumentDouble("-sa", 0.05, "Alpha for backtracking line search");
	step_beta_ = utils.ParseArgumentDouble("-sb", 0.3, "Beta for backtracking line search");
	MIN_F_ = 0.0;
	MAX_F_ = 1000.0;
	MIN_ETA_ = 0.0;
	MAX_ETA_ = 1000.0;
	MIN_P_ = 0.0001;
	MAX_P_ = 0.9999;
	cerr << "Model initiate done" << TimeString() << endl;
}

void Model::LoadGraph() {
	// system("dir");
	ifstream fin(input_graph_filename_);
	assert(fin && "open graph file failed");

	string s;
	int node_count, node, ref;
	fin >> node_count;
	getline(fin, s, '\n'); //remove the '\n' at the first line
	graph_.SetNodeCount(node_count);
	while (getline(fin, s, '\n')) {
		istringstream iss(s);
		iss >> node;
		map<int, int>::iterator it = id2order_.find(node);
		assert((it != id2order_.end()) && "No such paper in label file");
		node = it->second;
		while (iss >> ref) {
			it = id2order_.find(ref);
			assert((it != id2order_.end()) && "No such paper in label file");
			ref = it->second;
			graph_.AddEdge(ref, node);
		}
	}
	for (int i = 0; i < graph_.GetNodeCount(); ++i) {
		sort(graph_.GetInNeighbour(i).begin(), graph_.GetInNeighbour(i).end());
		sort(graph_.GetOutNeighbour(i).begin(), graph_.GetOutNeighbour(i).end());
	}
	cerr << "Load Graph Done. " << TimeString() << endl;
}

void Model::LoadLabel() {
	ifstream fin(input_label_filename_);
	assert(fin && "open label file failed");

	string s;
	int node_count, node;
	fin >> node_count;
	label_.resize(node_count);
	while (fin >> node) {
		label_[node].id = node;
		fin >> label_[node].date;
	}
	sort(label_.begin(), label_.end());
	for (int i = 0; i < label_.size(); ++i)
		id2order_.insert(pair<int, int>(label_[i].id, i));
	cerr << "Load Label Done. " << TimeString() << endl;
}


int Model::DetectCommunityCount() {
	return 100;
	//TODO: really detect optimal Community Number
}

int Model::GetOptimalCommunityCount() {
	if (optimal_community_count_ == -1)
		return DetectCommunityCount();
	return optimal_community_count_;
}

int Model::CountDifferent(vector<int> v1, vector<int> v2) {
	int count = 0, i1, i2;
	// In graph_init process, neighbours have been sorted. So it need not sort now.
	for (i1 = 0, i2 = 0; i1 < v1.size() && i2 < v2.size();)
		if (v1[i1] < v2[i2]) {
			++i1;
			++count;
		} else {
			if (v1[i1] == v2[i2])
				++i1;
			++i2;
		}
	count += v1.size() - i1;
	return count;
}

void Model::Resize() {
	F_.resize(graph_.GetNodeCount(), vector<double>(community_count_, 0));
	newF_.resize(graph_.GetNodeCount(), vector<double>(community_count_, 0));
	Eta_.resize(community_count_, vector<double>(community_count_, 1.0 / community_count_));
	newEta_.resize(community_count_, vector<double>(community_count_, 0));
	EtaF_.resize(graph_.GetNodeCount(), vector<double>(community_count_, 0));
	sumF_.insert(pair<int, vector<double>>(label_[0].date, vector<double>(community_count_, 0)));
	sumEtaF_.insert(pair<int, vector<double>>(label_[0].date, vector<double>(community_count_, 0)));
	backsumF_.insert(pair<int, vector<double>>(label_[0].date, vector<double>(community_count_, 0)));
	for (int i = 1; i < label_.size(); ++i) {
		assert((label_[i].date >= label_[i - 1].date) && "date label wrong");
		if (label_[i].date != label_[i - 1].date) {
			sumF_.insert(pair<int, vector<double>>(label_[i].date, vector<double>(community_count_, 0)));
			sumEtaF_.insert(pair<int, vector<double>>(label_[i].date, vector<double>(community_count_, 0)));
			backsumF_.insert(pair<int, vector<double>>(label_[i].date, vector<double>(community_count_, 0)));
		}
	}
}

void Model::InitNeighbourCommunity(int community_count) {
	community_count_ = community_count;
	Resize();

	//compute conductance of neighborhood community
	vector<pair<double, int>> conductances;
	conductances.reserve(graph_.GetNodeCount());
	for (int u = 0; u < graph_.GetNodeCount(); ++u) {
		// vector<int> &in_neighbour = graph_.GetInNeighbour(u);
		vector<int> S;
		S.push_back(u);
		S.insert(S.end(), graph_.GetInNeighbour(u).begin(), graph_.GetInNeighbour(u).end());
		// assert((S[0] < S[1]) && "ref error");
		int pos;
		for (pos = 1; pos < S.size() && S[pos] < u; ++pos)
			S[pos - 1] = S[pos];
		S[pos - 1] = u;
		int vol = 0, cut = 0;
		for (int i = 0; i < S.size(); ++i) {
			int v = S[i];
			vol += graph_.GetOutNeighbour(v).size();
			cut += CountDifferent(graph_.GetOutNeighbour(v), S);
		}
		conductances.push_back(pair<double, int>(vol != 0 ? ((double)cut / vol) : 10, u));
	}
	sort(conductances.begin(), conductances.end());
	cout << "conductance computation completed " << TimeString() << endl;
	ofstream fout("conductance.txt");
	for (int i = 0; i < conductances.size(); ++i)
		fout << conductances[i].first << '\t' << conductances[i].second << '\n';
	fout.close();

	//choose nodes with local minimum in conductance
	vector<bool> not_local_min(graph_.GetNodeCount(), false);
	int community_id = 0;
	for (int i = 0; i < conductances.size(); i++) {
		int u = conductances[i].second;
		if (not_local_min[u]) {
			continue;
		}
		//add the node and its neighbors to the current community
		F_[u][community_id] = 1;
		vector<int> &in_neighbour = graph_.GetInNeighbour(u);
		for (int i = 0; i < in_neighbour.size(); ++i)
			F_[in_neighbour[i]][community_id] = 1;
		//exclude its neighbors from the next considerations
		vector<int> &out_neighbour = graph_.GetOutNeighbour(u);
		for (int i = 0; i < out_neighbour.size(); ++i)
			not_local_min[out_neighbour[i]] = true;
		++community_id;
		if (community_id >= community_count) {
			break;
		}
	}
	if (community_count > community_id) {
		printf("%d communities needed to fill randomly\n", community_count - community_id);
	}
	//assign a member to zero-member community (if any)
	default_random_engine e(time(0));
	uniform_int_distribution<int> uniform_int(0, graph_.GetNodeCount() - 1);
	uniform_real_distribution<double> uniform_double(0, 1);
	for (; community_id < community_count; ++community_id)
		for (int i = 0; i < 10; ++i)
			F_[uniform_int(e)][community_id] += uniform_double(e);

	//TODO : Initiate eta.
}

void Model::CalculateEtaF(vector<vector<double>> &eta) {
	for (int u = 0; u < F_.size(); ++u)
		for (int i = 0; i < eta.size(); ++i)
			EtaF_[u][i] = VectorDot(eta[i], F_[u]);
}

void Model::AccumulateSumF() {
	vector<double> sum(community_count_, 0);
	int date = -1;
	map<int, vector<double>>::iterator it;
	for (int u = 0; u < F_.size(); ++u) {
		if (label_[u].date != date) {
			it = sumF_.find(label_[u].date);
			assert((it != sumF_.end()) && "No this date");
			it->second = sum;
			date = label_[u].date;
		}
		VectorAddTo(F_[u], sum);
	}
}

void Model::AccumulateSumEtaF() {
	vector<double> sum(community_count_, 0);
	int date = -1;
	map<int, vector<double>>::iterator it;
	for (int u = 0; u < EtaF_.size(); ++u) {
		if (label_[u].date != date) {
			it = sumEtaF_.find(label_[u].date);
			assert((it != sumEtaF_.end()) && "No this date");
			it->second = sum;
			date = label_[u].date;
		}
		VectorAddTo(EtaF_[u], sum);
	}
}

void Model::BackAccumulateSumF() {
	vector<double> sum(community_count_, 0);
	int date = -1;
	map<int, vector<double>>::iterator it;
	for (int u = F_.size() - 1; u >= 0; --u) {
		if (label_[u].date != date) {
			it = backsumF_.find(label_[u].date);
			assert((it != backsumF_.end()) && "No this date");
			it->second = sum;
			date = label_[u].date;
		}
		VectorAddTo(F_[u], sum);
	}
}

double Model::Norm2(vector<double> &v) {
	return VectorDot(v, v);
}

double Model::Norm2(vector<vector<double>> &m) {
	double result = 0.0;
	for (int i = 0; i < m.size(); ++i)
		result += VectorDot(m[i], m[i]);
	return result;
}

void Model::FGradientForRow(int u, vector<double> &gradient) {
	gradient.assign(community_count_, 0);
	for (int i = 0; i < graph_.GetInNeighbour(u).size(); ++i) {
		int v = graph_.GetInNeighbour(u)[i];
		double fu_eta_fv = VectorDot(F_[u], EtaF_[v]);
		double p_uv = min(max(exp(-fu_eta_fv), MIN_P_), MAX_P_);
		double argument = p_uv / (1 - p_uv) + 1; // this '1' means the second part of the fomulation
		for (int c = 0; c < community_count_; ++c)
			gradient[c] += argument * EtaF_[v][c];
		for (int i = 0; i < gradient.size(); ++i)
			if (isnan(gradient[i])) {
				cout << "nan occur!!" << argument << ' ' << fu_eta_fv;
			}
	}
	map<int, vector<double>>::iterator it = sumEtaF_.find(label_[u].date);
	assert((it != sumEtaF_.end()) && "No this date");
	VectorSubTo(it->second, gradient);
	for (int i = 0; i < gradient.size(); ++i)
		if (isnan(gradient[i])) {
			cout << "nan occur!!";
		}
}

double Model::FLikelihoodForRow(int u, vector<double> &Fu) {
	double L = 0.0;

	for (int i = 0; i < graph_.GetInNeighbour(u).size(); ++i) {
		int v = graph_.GetInNeighbour(u)[i];
		double fu_eta_fv = VectorDot(Fu, EtaF_[v]);
		double p_uv = min(max(exp(-fu_eta_fv), MIN_P_), MAX_P_);
		double Lu = log(1 - p_uv);
		L += Lu + fu_eta_fv;
	}
	map<int, vector<double>>::iterator it = sumEtaF_.find(label_[u].date);
	assert((it != sumEtaF_.end()) && "No this date");
	L -= VectorDot(Fu, it->second);

	vector<double> EtaFu(community_count_);
	for (int i = 0; i < community_count_; ++i)
		EtaFu[i] = VectorDot(Eta_[i], Fu);
	for (int i = 0; i < graph_.GetOutNeighbour(u).size(); ++i) {
		int v = graph_.GetOutNeighbour(u)[i];
		double fv_eta_fu = VectorDot(F_[v], EtaFu);
		double p_vu = min(max(exp(-fv_eta_fu), MIN_P_), MAX_P_);
		double Lu = log(1 - p_vu);
		L += Lu + fv_eta_fu;
	}
	it = backsumF_.find(label_[u].date);
	assert((it != backsumF_.end()) && "No this date");
	L -= VectorDot(it->second, EtaFu);
	return L;
}

bool Model::FGoStepSizeByLineSearch(int u, vector<double> &gradient, const int MAX_ITER) {
	double step_size = 1.0;
	double init_likelihood = FLikelihoodForRow(u, F_[u]);
	double square_gradient = VectorDot(gradient, gradient);
	for (int iter = 0; iter < MAX_ITER; iter++) {
		for (int c = 0; c < community_count_; c++) {
			double new_value = F_[u][c] + step_size * gradient[c];
			if (new_value < MIN_F_)
				new_value = MIN_F_;
			if (new_value > MAX_F_)
				new_value = MAX_F_;
			newF_[u][c] = new_value;
		}
		if (FLikelihoodForRow(u, newF_[u]) < init_likelihood + step_alpha_ * step_size * square_gradient)
			step_size *= step_beta_;
		else
			return true;
	}
	newF_[u] = F_[u];
	return 0;
}

void Model::UpdateF(vector<int> &order) {
	CalculateEtaF(Eta_);
	AccumulateSumEtaF();
	BackAccumulateSumF();
	vector<double> gradient;
	for (int i = 0; i < order.size(); ++i) {
		int u = order[i];
		FGradientForRow(u, gradient);
		//PrintVector("gradient:", gradient);
		//system("PAUSE");
		if (Norm2(gradient) < 1e-4)
			newF_[u] = F_[u];
		else
			FGoStepSizeByLineSearch(u, gradient);
	}
	F_ = newF_;
}

void Model::EtaGradient(vector<vector<double>> &gradient) {
	vector<double> tmp;
	for (int u = 0; u < graph_.GetNodeCount(); ++u) {
		tmp.assign(community_count_, 0);
		for (int i = 0; i < graph_.GetInNeighbour(u).size(); ++i) {
			int v = graph_.GetInNeighbour(u)[i];
			double fu_eta_fv = VectorDot(F_[u], EtaF_[v]);
			double p_uv = min(max(exp(-fu_eta_fv), MIN_P_), MAX_P_);
			double argument = p_uv / (1.00000000001 - p_uv) + 1; // this '1' means the second part of the fomulation
			for (int c = 0; c < community_count_; ++c)
				tmp[c] += argument * F_[v][c];
		}
		map<int, vector<double>>::iterator it = sumF_.find(label_[u].date);
		assert((it != sumF_.end()) && "No this date");
		VectorSubTo(it->second, tmp);
		for (int i = 0; i < community_count_; ++i)
			for (int j = 0; j < community_count_; ++j)
				gradient[i][j] = F_[u][i] * tmp[j];
	}
}

double Model::EtaLikelihood(vector<vector<double>> &eta) {
	CalculateEtaF(eta);
	AccumulateSumEtaF();

	double L = 0.0;
	for (int u = 0; u < graph_.GetNodeCount(); ++u) {
		for (int i = 0; i < graph_.GetInNeighbour(u).size(); ++i) {
			int v = graph_.GetInNeighbour(u)[i];
			double fu_eta_fv = VectorDot(F_[u], EtaF_[v]);
			double p_uv = min(max(exp(-fu_eta_fv), MIN_P_), MAX_P_);
			double Lu = log(1 - p_uv);
			L += Lu + fu_eta_fv;
		}
		map<int, vector<double>>::iterator it = sumEtaF_.find(label_[u].date);
		assert((it != sumEtaF_.end()) && "No this date");
		L -= VectorDot(F_[u], it->second);
	}
	return L;
}

bool Model::EtaGoStepSizeByLineSearch(vector<vector<double>> &gradient, const int MAX_ITER) {
	double step_size = 1.0;
	double init_likelihood = EtaLikelihood(Eta_);
	double square_gradient = 0.0;
	for (int c = 0; c < community_count_; ++c)
		square_gradient += VectorDot(gradient[c], gradient[c]);
	for (int iter = 0; iter < MAX_ITER; iter++) {
		for (int i = 0; i < community_count_; ++i)
			for (int j = 0; j < community_count_; ++j) {
				double new_value = Eta_[i][j] + step_size * gradient[i][j];
				if (new_value < MIN_ETA_)
					new_value = MIN_ETA_;
				if (new_value > MAX_ETA_)
					new_value = MAX_ETA_;
				newEta_[i][j] = new_value;
			}
		if (EtaLikelihood(newEta_) < init_likelihood + step_alpha_ * step_size * square_gradient)
			step_size *= step_beta_;
		else
			return true;
	}
	newEta_ = Eta_;
	return false;
}

void Model::UpdateEta() {
	CalculateEtaF(Eta_);
	AccumulateSumF();
	vector<vector<double>> gradient(community_count_, vector<double>(community_count_, 0));
	EtaGradient(gradient);
	if (Norm2(gradient) < 1e-4)
		return;
	EtaGoStepSizeByLineSearch(gradient);
	Eta_ = newEta_;
}

string Model::GenerateModelName(int iteration) {
	assert((iteration >= 0) && "iteration should be greater than zero");
	ostringstream ss;
	ss << setw(5) << iteration;
	string s = ss.str();
	replace(s.begin(), s.end(), ' ', '0');
	return s;
}

void Model::SaveF(string &model_name) {
	string filename = output_file_perfix_ + model_name + ".f";
	ofstream fout(filename);
	assert(fout && "open f ouput file failed");
	fout << F_.size() << '\t' << F_[0].size() << '\n';
	for (int i = 0; i < F_.size(); ++i) {
		for (int j = 0; j < F_[i].size(); ++j)
			fout << F_[i][j] << '\t';
		fout << '\n';
	}
	fout.close();
}

void Model::SaveEta(string &model_name) {
	string filename = output_file_perfix_ + model_name + ".eta";
	ofstream fout(filename);
	assert(fout && "open eta ouput file failed");
	fout << Eta_.size() << '\t' << Eta_[0].size() << '\n';
	for (int i = 0; i < Eta_.size(); ++i) {
		for (int j = 0; j < Eta_[i].size(); ++j)
			fout << Eta_[i][j] << '\t';
		fout << '\n';
	}
	fout.close();
}

void Model::SaveData(string model_name) {
	SaveF(model_name);
	SaveEta(model_name);
}

void Model::MLEGradAscent(const int MAX_ITER, const int RECORD_PERIOD) {
	int iter = 0;
	vector<int> order(graph_.GetNodeCount());
	for (int i = 0; i < order.size(); ++i)
		order[i] = i;
	default_random_engine e(time(0));
	for (iter = 0; iter < MAX_ITER; ++iter) {
		shuffle(order.begin(), order.end(), e);
		UpdateF(order);

		UpdateEta();

		if (iter % RECORD_PERIOD == 0) {
			cout << iter << " iterations. " << TimeString() << endl;
			SaveData(GenerateModelName(iter));
		}
	}
	SaveData("final");
}


// int main(int argc, char const *argv[]) {
// 	return 0;
// }
