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
	// vector<vector<int>> v;
	// v.resize(10, vector<int>(10));
	// for (int i = 0; i < 10; ++i)
	// 	for (int j = 0; j < 10; ++j)
	// 		v[i][j] = i * 10 + j;
	// for (int i = 0; i < 10; ++i) {
	// 	for (int j = 0; j < 10; ++j)
	// 		printf("%d ", v[i][j]);
	// 	printf("\n");
	// }
	// Test test;
	// printf("123\n");
	// const vector<Tmp> &v = test.Get(0);
	// // printf("%d\n", test.Get(0).size());
	// printf("321\n");
	// cout << TimeString();
	// default_random_engine e(time(0));
	// uniform_int_distribution<int> u(0,2);
	// for(int i = 0; i < 30; ++i)
	// 	printf("%d ", u(e));
	// ifstream fin("not_exist.txt");
	// cout << (bool)fin;
	// ifstream fin("test.txt");
	// string s;
	// while(getline(fin, s, '\n')){
	// 	printf("|%s|\n", s.c_str());
	// }
	// vector<vector<Tmp>> v1(3, vector<Tmp>(3)), v2;
	// vector<vector<Tmp>> v3(4, vector<Tmp>(3));
	// v3 = v1;
	// printf("------------\n");
	// v3.insert(v3.begin(), v1.begin(), v1.end());
	// v3.assign(v1.begin(), v1.end());
	// printf("%d %d\n", v3.size(), v3.capacity());
	// v3.resize(2);
	// int a[] = {1, 2, 3};
	// vector<int> v(a, a + 3);
	// v.clear();
	// printf("%d %d\n", v.size(), v.capacity());
	// v.resize(3);
	// for (auto i : v)
	// 	printf("%d ", i);
	string a = "sad";
	cout << a << flush;
	cout << '\r' << 'b';
}
