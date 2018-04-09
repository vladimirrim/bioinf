#include <iostream>
#include <fstream>
#include "../include/xmlParser.h"
#include "../include/massGeneration.h"
#include "../include/tsvParser.h"
#include <cmath>
#include <algorithm>

using namespace std;

const massTable MASS_TABLE;

void printAnnotatedPicks(const tsvParser &p, xmlParser &parser) {
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

                    if (abs(mass.first - prefix.second + MASS_TABLE.waterMass) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tprefix\twater loss";
                    }

                    if (abs(mass.first - prefix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
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

                    if (abs(mass.first - prefix.second + MASS_TABLE.waterMass) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix\twater loss";
                    }

                    if (abs(mass.first - prefix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
                        if (!isCovered[prefix.first]) {
                            covered++;
                            isCovered[prefix.first] = true;
                            isCovered[ts.linkedPairs[prefix.first]] = true;
                        }
                        os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix\tammonia loss";
                    }
                }
            }

            os << "\n covered " << (double) covered / (ts.sequenceLength - 1) * 100 << "% of peptide links\n\n";
        }
    }
    os.close();
}

void showAnnotatedPicks(ofstream &os, TheorySpectra &ts, MSResult &spectra, xmlParser &parser) {
    os << "\n" << spectra.scanId << "\t" << spectra.eValue << "\t" << (spectra.aType == HCD ? "HCD" : "CID");
    int covered = 0;
    map<string, bool> isCovered;
    double epsilon = 0.02;
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

            if (abs(mass.first - prefix.second + MASS_TABLE.waterMass) <= epsilon) {
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tprefix\twater loss";
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
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

            if (abs(mass.first - prefix.second + MASS_TABLE.waterMass) <= epsilon) {
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix\twater loss";
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix\tammonia loss";
            }
        }
    }

    os << "\n covered " << (double) covered / (ts.sequenceLength - 1) * 100 << "% of peptide links\n\n";
}

double calculateAnnotatedPicks(const MSResult &spectra, const vector<MSResult> &elem, TheorySpectra &ts,
                               xmlParser &parser) {
    double epsilon = 0.02;
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
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.waterMass) <= epsilon) {
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
            }

        }
    }
    for (auto &suffix:ts.suffixes) {
        string name = suffix.first;
        for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
            if (abs(mass.first - suffix.second) <= epsilon) {
                if (!isCovered[suffix.first]) {
                    covered++;
                    isCovered[suffix.first] = true;
                    isCovered[ts.linkedPairs[suffix.first]] = true;
                }
            }

            if (abs(mass.first - suffix.second + MASS_TABLE.waterMass) <= epsilon) {
                if (!isCovered[suffix.first]) {
                    covered++;
                    isCovered[suffix.first] = true;
                    isCovered[ts.linkedPairs[suffix.first]] = true;
                }
            }

            if (abs(mass.first - suffix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
                if (!isCovered[suffix.first]) {
                    covered++;
                    isCovered[suffix.first] = true;
                    isCovered[ts.linkedPairs[suffix.first]] = true;
                }
            }
        }
    }
    return (double) covered / (ts.sequenceLength - 1) * 100;
}

void printNeutralAndPlusHDifference(const tsvParser &p, xmlParser &parser, const massTable &t) {
    ofstream os("difference.txt");
    vector<pair<double, MSResult>> differences;
    for (auto &elem: p.results) {
        string sequence = elem.first;
        TheorySpectra tsNeutral = generateMasses(sequence);
        TheorySpectra tsPlusH = generateMasses(sequence, t.HMass);
        for (auto &spectra: elem.second) {
            double coveredNeutral = calculateAnnotatedPicks(spectra, elem.second, tsNeutral, parser);
            double coveredPlusH = calculateAnnotatedPicks(spectra, elem.second, tsPlusH, parser);
            differences.emplace_back(abs(coveredNeutral - coveredPlusH), spectra);
        }
    }
    sort(differences.begin(), differences.end(),
         [](pair<double, MSResult> x, pair<double, MSResult> y) { return x.first < y.first; });
    for (int i = 1; i <= 3; i++) {
        MSResult spectra = differences[differences.size() - i].second;
        string sequence = spectra.sequence;
        TheorySpectra tsNeutral = generateMasses(sequence);
        TheorySpectra tsPlusH = generateMasses(sequence, t.HMass);
        os << sequence << "\n";
        os << "difference: " << differences[differences.size() - i].first << "\n";
        os << "Neutral:\n";
        showAnnotatedPicks(os, tsNeutral, spectra, parser);
        os << "Plus H: \n";
        showAnnotatedPicks(os, tsPlusH, spectra, parser);
    }
    os.close();
}

int main() {
    xmlParser parser;
    parser.parseXML("../120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.mzXML.good.msalign");
    tsvParser p;
    p.parseTSV("../120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.tsv");
    printNeutralAndPlusHDifference(p, parser, MASS_TABLE);
    return 0;
}