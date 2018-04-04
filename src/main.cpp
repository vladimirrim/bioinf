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
    double epsilon = 0.02;
    for (auto &elem: p.results) {
        string sequence = elem.first;
        os << sequence << "\n";
        TheorySpectra ts = generateMasses(sequence);
        for (auto &spectra: elem.second) {
            int cnt = 0;
            os << "\n" << spectra.scanId << "\t" << spectra.eValue << "\t" << (spectra.aType == HCD ? "HCD" : "CID");
            long double sum = 0;
            long double covered = 0;
            for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities)
                sum += mass.second;
            for (auto &prefix:ts.prefixes) {
                string name = prefix.first;
                for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
                    if (abs(mass.first - prefix.second) <= epsilon) {
                        cnt++;
                        covered += mass.second;
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tprefix";
                    }
                }
            }
            for (auto &prefix:ts.suffixes) {
                string name = prefix.first;
                for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
                    if (abs(mass.first - prefix.second) <= epsilon) {
                        cnt++;
                        covered += mass.second;
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix";
                    }
                }
            }
            os << "\n cover " << (double) cnt / (parser.spectras[spectra.scanId].massesAndIntensities.size()) * 100
               << " percent -- by picks count";
            long double res = covered / sum * 100;
            os << "\n cover " << (double) res << " percent -- by intensity";
            os << "\n";
        }
    }
    os.close();
    return 0;
}