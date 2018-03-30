#include <iostream>
#include <fstream>
#include "../include/xmlParser.h"
#include "../include/massGeneration.h"
#include "../include/tsvParser.h"
#include <cmath>

using namespace std;

int main() {
    xmlParser parser;
    parser.parseXML("../120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.good.msalign");
    massTable t;
    tsvParser p;
    p.parseTSV("../120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.tsv");
    ofstream os("results.txt");
    double epsilon = 3;
    for (auto &elem: p.results) {
        string sequence = elem.first;
        os << sequence << "\n";
        TheorySpectra ts = generateMasses(sequence);
        for (auto &spectra: elem.second) {
            int cnt =0;
            os << "\n" << spectra.scanId << "\t" << spectra.eValue << "\t" << (spectra.aType == HCD ? "HCD" : "CID");
            for (auto &prefix:ts.prefixes) {
                string name = prefix.first;
                for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
                    if (abs(mass.first - prefix.second) <= epsilon) {
                        cnt++;
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tprefix";
                    }
                }
            }
            for (auto &prefix:ts.suffixes) {
                string name = prefix.first;
                for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
                    if (abs(mass.first - prefix.second) <= epsilon) {
                        cnt++;
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix";
                    }
                }
            }
            os << "\n cover " << (double) cnt/(parser.spectras[spectra.scanId].massesAndIntensities.size()) << " percent";
            os << "\n";
        }
    }
    os.close();
    return 0;
}