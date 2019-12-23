#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

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
	for (unsigned int i = 0; targetSize <= saRNA.size()-i; ++i) {
		std::string target = saRNA.substr(i, targetSize);
		
		if (!gcContent(target)) continue;
		if (!consecutive(target)) continue;
		if (!isNucleotideXY(target, 0, 'G', 'C')) continue;
		if (!isNucleotideXY(target, 1, 'G', 'C')) continue;
		if (!isNucleotideXY(target, targetSize-2, 'A', 'T')) continue;
		if (!isNucleotideXY(target, targetSize-1, 'A', 'A')) continue;

		if (saRNA.size()-i >= targetSize+2) {
			if (!isNucleotideXY(saRNA, i+targetSize+1, 'A', 'T')) continue;
			if (!isNucleotideXY(saRNA, i+targetSize, 'A', 'T')) continue;
		} else if (saRNA.size()-i == targetSize+1) {
			if (!isNucleotideXY(saRNA, i+targetSize, 'A', 'T')) continue;
		}

		std::cout << target << '\n';
	}
}

// usage(argv[]) displays usage error message and terminates program.
void usage(char *argv[]) {
	std::cerr << "Usage: " << argv[0] << " [ input-file ]" << std::endl;
	std::exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
	switch ( argc ) {
		case 2:
			break;
		default:
			usage(argv);
	}

	std::ifstream in( argv[1] );
	std::string saRNA;

	for ( ;; ) {
		in >> saRNA;
		if ( in.fail() ) break;
	}

	rnaIter(saRNA, 19);
}
