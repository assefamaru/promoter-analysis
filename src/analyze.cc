#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>
#include <unordered_map>

// gcContent(target) returns true if the target GC Content
// is between 40-60%, and false otherwise.
bool gcContent(std::string target) {
	size_t numG = std::count(target.begin(), target.end(), 'G');
	size_t numC = std::count(target.begin(), target.end(), 'C');
	double gc = ((numG + numC) / 19.0) * 100.0;

	return 40.0 <= gc && gc <= 60.0;
}

// consecutive(target) returns true if the target contains
// (<= 3) consecutive nucleotides, and false otherwise.
bool consecutive(std::string target) {
	return target.find("AAAA") == std::string::npos &&
		target.find("TTTT") == std::string::npos &&
		target.find("GGGG") == std::string::npos &&
		target.find("CCCC") == std::string::npos;
}

// triRepeats(target) returns the number of tri-repeats
// (or 3 consecutive nucleotides) in the target.
int triRepeats(std::string target) {
	int count = 0;

	std::string::size_type pos = 0;
	std::string sub = "AAA";
	while ((pos = target.find(sub, pos)) != std::string::npos) {
		++count;
		pos += sub.size();
	}

	pos = 0;
	sub = "TTT";
	while ((pos = target.find(sub, pos)) != std::string::npos) {
		++count;
		pos += sub.size();
	}

	pos = 0;
	sub = "GGG";
	while ((pos = target.find(sub, pos)) != std::string::npos) {
		++count;
		pos += sub.size();
	}

	pos = 0;
	sub = "CCC";
	while ((pos = target.find(sub, pos)) != std::string::npos) {
		++count;
		pos += sub.size();
	}

	return count;
}

// deltaG(seq, len) returns the delta G calculated using the nearest-neighbor method.
// Table and method can be found at: 
// https://en.wikipedia.org/wiki/Nucleic_acid_thermodynamics#Nearest-neighbor_method
double deltaG(std::string seq, int len) {
	double deltaGTotal = 0.0;

	std::unordered_map<std::string, double> umap;

	umap["AA"] = -4.26;
	umap["TT"] = -4.26;
	umap["AT"] = -3.67;
	umap["TA"] = -2.50;
	umap["CA"] = -6.12;
	umap["TG"] = -6.12;
	umap["GT"] = -6.09;
	umap["AC"] = -6.09;
	umap["AG"] = -5.40;
	umap["CT"] = -5.40;
	umap["GA"] = -5.51;
	umap["TC"] = -5.51;
	umap["CG"] = -9.07;
	umap["GC"] = -9.36;
	umap["GG"] = -7.66;
	umap["CC"] = -7.66;
	umap["A"] = 4.31;
	umap["T"] = 4.31;
	umap["G"] = 4.05;
	umap["C"] = 4.05;

	deltaGTotal += umap[seq.substr(0, 1)];
	for (int i = 0; i < len-1; ++i) {
		deltaGTotal += umap[seq.substr(i, 2)];
	}
	deltaGTotal += umap[seq.substr(len-1, 1)];

	return deltaGTotal;
}

// isNucleotideXY(target, index, X, Y) returns true if the target saRNA
// at the specified index contains either X or Y nucleotides, and false
// otherwise.
bool isNucleotideXY(std::string target, int index, char X, char Y) {
	return target[index] == X || target[index] == Y;
}

// rnaIter(saRNA, targetSize) iterates over the saRNA string,
// creating target substrings of size == targetSize.
// Each target saRNA is analyzed using various criteria.
// The criteria are as follows:
// * If target's GC Content is not between 40%-60%, it's filtered out.
// * If target contains > 3 consecutive nucleotides, it's filtered out.
// * If the ΔG (delta G) of target's left end is >= right end's, it's filtered out.
// - If 1st nucleotide of target is 'G' or 'C', its rank is increased by 10.
// - If 2nd nucleotide of target is 'G' or 'C', its rank is increased by 10.
// - If 18th nucleotide of target is 'A' or 'T', its rank is increased by 10.
// - If 19th nucleotide of target is 'A', its rank is increased by 10.
// - If 19th nucleotide of target is 'T', its rank is increased by 9.
// - If target contains tri-repeats (3 consecutive nucleotides), then its
//   rank is reduced as: rank = rank - ((# of tri-repeats) * 10)
// - If 20th nucleotide is 'A' or 'T', target's rank is increased by 4.
// - If 21st nucleotide is 'A' or 'T', target's rank is increased by 3.
// - If 22nd nucleotide is 'A' or 'T', target's rank is increased by 2.
// - If 23rd nucleotide is 'A' or 'T', target's rank is increased by 1.
// After filtering out and ranking 19nt targets, each target is printed to stdout,
// along with the 20-23 nucleotides, and their overall rank in descending order. 
void rnaIter(std::string saRNA, unsigned int targetSize) {
	struct target {
		int rank;
		std::string sequence;
		std::string outer;

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
		if (deltaG(current.substr(0, 4), 4) >= deltaG(current.substr(14, 4), 4)) continue;

		targets.push_back(target());
		targets[j].rank = 0;
		targets[j].sequence = current;
		targets[j].outer = "";

		if (isNucleotideXY(current, 0, 'G', 'C')) targets[j].rank += 10;
		if (isNucleotideXY(current, 1, 'G', 'C')) targets[j].rank += 10;
		if (isNucleotideXY(current, targetSize-2, 'A', 'T')) targets[j].rank += 10;
		if (isNucleotideXY(current, targetSize-1, 'A', 'A')) targets[j].rank += 10;
		if (isNucleotideXY(current, targetSize-1, 'T', 'T')) targets[j].rank += 9;
		targets[j].rank -= (triRepeats(current) * 10);

		for (int k = 4; k > 0; --k) {
			if (saRNA.size()-i >= targetSize+k) {
				if (targets[j].outer == "") {
					targets[j].outer = saRNA.substr(i+targetSize, k);
				}
				if (isNucleotideXY(saRNA, i+targetSize+(k-1), 'A', 'T')) {
					targets[j].rank += (4-(k-1));
				}
			}
		}

		targets[j].outer.resize (4, '-');

		++j;
	}

	std::sort(targets.begin(), targets.end(), std::greater<target>());

	for (target x : targets) {
		std::cout << x.sequence << "   " << x.outer << "   " << x.rank << '\n';
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
