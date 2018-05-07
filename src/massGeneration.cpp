#include <vector>
#include "../include/massGeneration.h"


using namespace std;

TheorySpectra generateMasses(const string &peptide, double additionalMass) {
    string prefix;
    massTable table;
    long double mass = additionalMass;
    int sequenceLength = 1;
    vector<pair<string, long double>> prefixes;
    map<string, string> linkedPairs;
    for (int i = 0; i < peptide.length() - 1; i++) {
        sequenceLength++;
        if (peptide[i] == 'C')
            mass += 57.021;
        mass += table[peptide[i]];
        prefix += peptide[i];
        prefixes.emplace_back(prefix, mass);
        if (peptide[i] == 'C')
            i += 7;
    }
    double globalMass = mass + table[peptide[peptide.length() - 1]];
    string suffix;
    mass = table.waterMass + additionalMass;
    vector<pair<string, long double>> suffixes;
    for (int i = peptide.length() - 1; i > 0; i--) {
        if (peptide[i] == '1') {
            mass += 57.021;
            mass += table['C'];
            suffix.insert(0, 1, 'C');
            i -= 7;
        } else {
            mass += table[peptide[i]];
            suffix.insert(0, 1, peptide[i]);
        }
        suffixes.emplace_back(suffix, mass);
    }
    for (int i = 0; i < sequenceLength - 1; i++) {
        linkedPairs[suffixes[i].first] = prefixes[sequenceLength - i - 2].first;
        linkedPairs[prefixes[sequenceLength - i - 2].first] = suffixes[i].first;
    }
    return TheorySpectra(sequenceLength, peptide, prefixes, suffixes, linkedPairs, globalMass);
}

double calculatePeptideMass(std::string peptide) {
    massTable table;
    double mass = 0;
    for (char i : peptide) {
        if (i == 'C')
            mass += 57.021;
        mass += table[i];
    }
    return mass;
}
