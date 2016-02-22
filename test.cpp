#include <string>
#include <cstdio>
#include <exception>
#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

struct Tmp {
	Tmp() {
		printf("construct\n");
	}
	Tmp(const Tmp &tmp2) {
		printf("copy construct\n");
	}
	Tmp &operator=(const Tmp &tmp2) {
		printf("=\n");
	}
	~Tmp() {
		printf("destruct\n");
	}
};

class Test {
public:
	Test() {
		v.resize(2, vector<Tmp>(2));
	}
	vector<Tmp> &Get(int num) {
		return v[num];
	}
private:
	vector<vector<Tmp>> v;
};


string TimeString() {
	time_t timep;
	time(&timep);
	return ctime(&timep);
}

int CountDifferent(vector<int> v1, vector<int> v2) {
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

int main(int argc, char const *argv[]) {
	int n = 10;
	char* buffer = new char[n];
	buffer[n] = 0;
}
