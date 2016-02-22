#ifndef COMMON_H
#define COMMON_H
#include <string>
#include <vector>
#include <ctime>
using namespace std;
inline string TimeString() {
	time_t timep;
	time(&timep);
	string t = ctime(&timep);
	return t.substr(0, t.length() - 1);
}

template<typename T>
T VectorDot(const vector<T> &v1, const vector<T> &v2) {
	assert((v1.size() == v2.size()) && "vector's length not match");
	T result = 0;
	for (int i = 0; i < v1.size(); ++i)
		result += v1[i] * v2[i];
	return result;
}

template<typename T>
vector<T> VectorAdd(const vector<T> &v1, const vector<T> &v2) {
	assert((v1.size() == v2.size()) && "vector's length not match");
	vector<T> result(v1.size());
	for (int i = 0; i < v1.size(); ++i)
		result[i] = v1[i] + v2[i];
	return result;
}

template<typename T>
vector<T> VectorAdd(const vector<T> &v1, const vector<T> &v2, vector<T> &v3) {
	assert((v1.size() == v2.size()) && "vector's length not match");
	for (int i = 0; i < v1.size(); ++i)
		v3[i] = v1[i] + v2[i];
}

template<typename T>
void VectorAddTo(const vector<T> &v1, vector<T> &v2) {
	assert((v1.size() == v2.size()) && "vector's length not match");
	for (int i = 0; i < v1.size(); ++i)
		v2[i] += v1[i];
}

template<typename T>
void VectorSubTo(const vector<T> &v1, vector<T> &v2) {
	assert((v1.size() == v2.size()) && "vector's length not match");
	for (int i = 0; i < v1.size(); ++i)
		v2[i] -= v1[i];
}

#endif