#include <vector>
#include "../include/massGeneration.h"


using namespace std;

TheorySpectra generateMasses(const string &peptide) {
    string prefix;
    massTable table;
    long double mass = 1.007;
    vector<pair<string, long double>> prefixes;
    for (int i = 0; i < peptide.length() - 1; i++) {
        if (peptide[i] == 'C')
            mass += 57.021;
        mass += table[peptide[i]];
        prefix += peptide[i];
        prefixes.emplace_back(prefix, mass);
        if (peptide[i] == 'C')
            i += 7;
    }
    string suffix;
    mass = table.waterMass + 1.007;
    vector<pair<string, long double>> suffixes;
    for (int i = peptide.length() - 1; i > 0; i--) {
        if (peptide[i] == '1') {
            mass += 57.021;
            mass += table['C'];
            suffix.insert(0, 1, 'C');
            i -= 7;
        }else{
            mass += table[peptide[i]];
            suffix.insert(0, 1, peptide[i]);
        }
        suffixes.emplace_back(suffix, mass);
    }
    return TheorySpectra(peptide, prefixes, suffixes);
}