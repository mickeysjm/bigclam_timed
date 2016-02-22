#include <string>
#include <vector>
#include <set>
#include <map>
#include <exception>
#include "utils.h"
#include "model.h"
using namespace std;

int main(int argc, char const *argv[]) {
	try {
		Utils utils(argc, argv);
		Model model;
		model.Init(utils);
		model.LoadGraph();
		model.LoadLabel();
		model.InitNeighbourCommunity(model.GetOptimalCommunityCount());
		model.SaveData("init");
		model.MLEGradAscent();
	} catch (exception &e) {
		cerr << "Exception caught: " << e.what() << endl;
	}
	system("PAUSE");
	return 0;
}