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
            os << "\n" << spectra.scanId << "\t" << spectra.eValue << "\t" << (spectra.aType == HCD ? "HCD" : "CID");
            int covered = 0;
            map<string, bool> isCovered;
            for (auto &prefix:ts.prefixes) {
                string name = prefix.first;
                for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
                    if (abs(mass.first - prefix.second) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tprefix";
                    }

                    if (abs(mass.first - prefix.second + t.waterMass) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tprefix\twater loss";
                    }

                    if (abs(mass.first - prefix.second + t.ammoniaMass) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tprefix\tammonia loss";
                    }

                }
            }
            for (auto &prefix:ts.suffixes) {
                string name = prefix.first;
                for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
                    if (abs(mass.first - prefix.second) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix";
                    }

                    if (abs(mass.first - prefix.second + t.waterMass) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix\twater loss";
                    }

                    if (abs(mass.first - prefix.second + t.ammoniaMass) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix\tammonia loss";
                    }
                }
            }

            os << "\n covered " << (double) covered / (ts.sequence.length() - 1) * 100 << "% of peptide links\n\n";
        }
    }
    os.close();
    return 0;
}