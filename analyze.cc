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

// nucl1GorC(target) returns true if the 1st nucleotide of the 
// target saRNA is either G or C, and false otherwise.
bool nucl1GorC(std::string target) {
	return target[0] == 'G' || target[0] == 'C';
}

// nucl2GorC(target) returns true if the 2nd nucleotide of the
// target saRNA is either G or C, and false otherwise.
bool nucl2GorC(std::string target) {
	return target[1] == 'G' || target[1] == 'C';
}

// nucl18AorT(target) returns true if the 18th nucleotide of the
// target saRNA is either A or T, and false otherwise.
bool nucl18AorT(std::string target) {
	return target[17] == 'A' || target[17] == 'T';
}

// nucl19A(target) returns true if the 19th nucleotide of the 
// target saRNA is A, and false otherwise.
bool nucl19A(std::string target) {
	return target[18] == 'A';
}

// rnaIter(saRNA, targetSize) iterates over the saRNA string,
// creating target substrings of size == targetSize. Each
// target rna is analyzed using various criteria of different 
// priority, in order to filter ideal targets in the saRNA.
void rnaIter(std::string saRNA, unsigned int targetSize) {
	for (unsigned int i = 0; targetSize < saRNA.size()-i; ++i) {
		std::string target = saRNA.substr(i, targetSize);

		if (gcContent(target) && consecutive(target) && nucl1GorC(target) && nucl2GorC(target) && nucl18AorT(target) && nucl19A(target)) {
			std::cout << target << '\n';
		}
	}
}

// usage(argv[]) displays usage error message and terminates program.
void usage(char *argv[]) {
	std::cerr << "Usage: " << argv[0] << " [ input-file ]" << std::endl;
	std::exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
	std::ifstream in( argv[1] );
	std::string saRNA;
	
	switch (argc) {
		case 2:
			for ( ;; ) {
				in >> saRNA;
				if ( in.fail() ) break;
			}
			break;
		default:
			usage(argv);
	}

	rnaIter(saRNA, 19);
}
