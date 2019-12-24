#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>

// gcContent(target) returns true if the target GC Content
// is between 40-60%, and false otherwise.
bool gcContent(std::string target) {
	size_t numG = std::count(target.begin(), target.end(), 'G');
	size_t numC = std::count(target.begin(), target.end(), 'C');
	double gc = ((numG + numC) / 19.0) * 100.0;

	return 40.0 <= gc && gc <= 60.0;
}

// consecutive(target) returns true if the target contains
// (<= 5) consecutive nucleotides, and false otherwise.
bool consecutive(std::string target) {
	return target.find("AAAAAA") == std::string::npos &&
		target.find("TTTTTT") == std::string::npos &&
		target.find("GGGGGG") == std::string::npos &&
		target.find("CCCCCC") == std::string::npos;
}

// isNucleotideXY(target, index, X, Y) returns true if the target saRNA
// at the specified index contains either X or Y nucleotides, and false
// otherwise.
bool isNucleotideXY(std::string target, int index, char X, char Y) {
	return target[index] == X || target[index] == Y;
}

// rnaIter(saRNA, targetSize) iterates over the saRNA string,
// creating target substrings of size == targetSize. Each
// target rna is analyzed using various criteria of different 
// priority, in order to filter ideal targets in the saRNA.
void rnaIter(std::string saRNA, unsigned int targetSize) {
	struct target {
		int rank;
		std::string sequence;

		bool operator>(const target &a) const {
			return rank > a.rank;
		}
	};

	std::vector<target> targets;
	
	int j = 0;
	for (unsigned int i = 0; targetSize <= saRNA.size()-i; ++i) {
		std::string current = saRNA.substr(i, targetSize);

		if (!gcContent(current)) continue;
		if (!consecutive(current)) continue;
		
		targets.push_back(target());
		targets[j].rank = 0;
		targets[j].sequence = current;

		if (isNucleotideXY(current, 0, 'G', 'C')) ++targets[j].rank;
		if (isNucleotideXY(current, 1, 'G', 'C')) ++targets[j].rank;
		if (isNucleotideXY(current, targetSize-2, 'A', 'T')) ++targets[j].rank;
		if (isNucleotideXY(current, targetSize-1, 'A', 'A')) ++targets[j].rank;

		if (saRNA.size()-i >= targetSize+2) {
			if (isNucleotideXY(saRNA, i+targetSize+1, 'A', 'T')) ++targets[j].rank;
			if (isNucleotideXY(saRNA, i+targetSize, 'A', 'T')) ++targets[j].rank;
		} else if (saRNA.size()-i == targetSize+1) {
			if (isNucleotideXY(saRNA, i+targetSize, 'A', 'T')) ++targets[j].rank;
		}

		++j;
	}

	std::sort(targets.begin(), targets.end(), std::greater<target>());
	
	for (target x : targets) {
		std::cout << x.sequence << '\n';
	}
}

// error(message) displays error message and terminates program.
void error(std::string message) {
	std::cerr << message << std::endl;
	std::exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
	switch ( argc ) {
		case 2:
			break;
		default:
			std::string name = argv[0];
			error("Usage: " + name + " [ input-file ]");
	}

	std::ifstream in( argv[1] );
	std::string saRNA;

	for ( ;; ) {
		in >> saRNA;
		if ( in.fail() ) break;
	}

	if (saRNA.size() < 19) {
		error("saRNA must contain at least 19 nucleotides!");
	} else {
		rnaIter(saRNA, 19);
	}
}
