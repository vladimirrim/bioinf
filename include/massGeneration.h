#pragma once


#include <map>
#include <utility>

class massTable {
public:
    std::map<char, long double> masses;
    const long double waterMass = 18.010565;

    massTable() {
        masses['A'] = 71.03711;
        masses['R'] = 156.10111;
        masses['N'] = 114.04293;
        masses['D'] = 115.02694;
        masses['C'] = 103.00919;
        masses['E'] = 129.04259;
        masses['Q'] = 128.05858;
        masses['G'] = 57.02146;
        masses['H'] = 137.05891;
        masses['I'] = 113.08406;
        masses['L'] = 113.08406;
        masses['K'] = 128.09496;
        masses['M'] = 131.04049;
        masses['F'] = 147.06841;
        masses['P'] = 97.05276;
        masses['S'] = 87.03203;
        masses['T'] = 101.04768;
        masses['W'] = 186.07931;
        masses['Y'] = 163.06333;
        masses['V'] = 99.06841;

    }

    long double operator[](char c) {
        return masses[c];
    }
};

struct TheorySpectra {
    std::string sequence;
    std::vector<std::pair<std::string, long double>> prefixes{};
    std::vector<std::pair<std::string, long double>> suffixes{};

    TheorySpectra(std::string sequence, std::vector<std::pair<std::string, long double>> prefixes,
                  std::vector<std::pair<std::string, long double>> suffixes) :
            sequence(std::move(sequence)), prefixes(prefixes), suffixes(suffixes) {};
};

TheorySpectra generateMasses(const std::string &peptide);