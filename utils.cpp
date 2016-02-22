#include "utils.h"
#include "common.h"
#include <cassert>

Utils::Utils(int argc, char const *argv[]) {
	assert((argc & 1) && "argument count error");
	for (int i = 1; i < argc; i += 2) {
		arguments_.insert(pair<string, string>(argv[i], argv[i + 1]));
	}
	cerr << "Read arguments done" << TimeString() << endl;
	cout << "usage:bigclam.exe" << endl;
}

int Utils::ParseArgumentInt(string identifier, int init_argument, string usage){
	cout << '\t' << identifier << ':' << usage << "(default" << init_argument << ")\n";
	map<string, string>::iterator it = arguments_.find(identifier);
	if(it == arguments_.end())
		return init_argument;
	return stoi(it->second);
}

string Utils::ParseArgumentString(string identifier, string init_argument, string usage){
	cout << '\t' << identifier << ':' << usage << "(default" << init_argument << ")\n";
	map<string, string>::iterator it = arguments_.find(identifier);
	if(it == arguments_.end())
		return init_argument;
	return it->second;	
}

double Utils::ParseArgumentDouble(string identifier, double init_argument, string usage){
	cout << '\t' << identifier << ':' << usage << "(default" << init_argument << ")\n";
	map<string, string>::iterator it = arguments_.find(identifier);
	if(it == arguments_.end())
		return init_argument;
	return stod(it->second);
}

// int main(int argc, char const *argv[])
// {
	
// 	return 0;
// }