
#include "sr_apx/sr_apx.hpp"

int main(int argc, char* argv[]) {
	std::string filename = argv[1];
	int index = filename.find(".txt");
	std::string output = filename.substr(0, index) + ".s6";

	sr_apx::Graph g = sr_apx::read_edge_list(filename.c_str());
	sr_apx::Graph h = sr_apx::shuffle_vertices(g);
	sr_apx::write_sparse6(g, output.c_str());
}
