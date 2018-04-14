#include <iostream>
#include <fstream>
#include "../include/xmlParser.h"
#include "../include/massGeneration.h"
#include "../include/tsvParser.h"
#include <cmath>
#include <algorithm>

using namespace std;

const massTable MASS_TABLE;

pair<double, double> calculateError(double theoryMass, double mass, double globalMass) {
    double absError = abs(mass - theoryMass);
    double relError = absError / globalMass * 100;
    return {absError, relError};
}

void showAnnotatedPicks(ofstream &os, TheorySpectra &ts, const MSResult &spectra, xmlParser &parser) {
    os << "\n" << spectra.scanId << "\t" << spectra.eValue << "\t" << (spectra.aType == HCD ? "HCD" : "CID");
    int covered = 0;
    map<string, bool> isCovered;
    double epsilon = 2;
    int annotatedCnt = 0;
    double averageRelError = 0, averageAbsError = 0, maxRelError = 0, maxAbsError = 0, minRelError = 1, minAbsError = 1;
    for (auto &prefix:ts.prefixes) {
        string name = prefix.first;
        for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
            if (abs(mass.first - prefix.second) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                pair<double, double> errors = calculateError(prefix.second, mass.first, ts.globalMass);
                averageAbsError += errors.first;
                averageRelError += errors.second;
                maxAbsError = max(maxAbsError, errors.first);
                maxRelError = max(maxRelError, errors.second);
                minAbsError = min(minAbsError, errors.first);
                minRelError = min(minRelError, errors.second);
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tprefix\t\tabsolute mass error: "
                   << errors.first <<
                   "\trelative mass error: " << errors.second << "%";
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.waterMass) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                pair<double, double> errors = calculateError(prefix.second - MASS_TABLE.waterMass, mass.first,
                                                             ts.globalMass - MASS_TABLE.waterMass);
                averageAbsError += errors.first;
                averageRelError += errors.second;
                maxAbsError = max(maxAbsError, errors.first);
                maxRelError = max(maxRelError, errors.second);
                minAbsError = min(minAbsError, errors.first);
                minRelError = min(minRelError, errors.second);
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                   << "\tprefix\twater loss\tabsolute mass error: " << errors.first <<
                   "\trelative mass error: " << errors.second << "%";
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                pair<double, double> errors = calculateError(prefix.second - MASS_TABLE.ammoniaMass, mass.first,
                                                             ts.globalMass - MASS_TABLE.ammoniaMass);
                averageAbsError += errors.first;
                averageRelError += errors.second;
                maxAbsError = max(maxAbsError, errors.first);
                maxRelError = max(maxRelError, errors.second);
                minAbsError = min(minAbsError, errors.first);
                minRelError = min(minRelError, errors.second);
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                   << "\tprefix\tammonia loss\tabsolute mass error: " << errors.first <<
                   "\trelative mass error: " << errors.second << "%";
            }

        }
    }
    for (auto &prefix:ts.suffixes) {
        string name = prefix.first;
        for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
            if (abs(mass.first - prefix.second) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                pair<double, double> errors = calculateError(prefix.second, mass.first, ts.globalMass);
                averageAbsError += errors.first;
                averageRelError += errors.second;
                maxAbsError = max(maxAbsError, errors.first);
                maxRelError = max(maxRelError, errors.second);
                minAbsError = min(minAbsError, errors.first);
                minRelError = min(minRelError, errors.second);
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name << "\tsuffix\t\tabsolute mass error: "
                   << errors.first <<
                   "\trelative mass error: " << errors.second << "%";
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.waterMass) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                pair<double, double> errors = calculateError(prefix.second - MASS_TABLE.waterMass, mass.first,
                                                             ts.globalMass - MASS_TABLE.waterMass);
                averageAbsError += errors.first;
                averageRelError += errors.second;
                maxAbsError = max(maxAbsError, errors.first);
                maxRelError = max(maxRelError, errors.second);
                minAbsError = min(minAbsError, errors.first);
                minRelError = min(minRelError, errors.second);
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                   << "\tsuffix\twater loss\tabsolute mass error: " << errors.first <<
                   "\trelative mass error: " << errors.second << "%";
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                pair<double, double> errors = calculateError(prefix.second - MASS_TABLE.ammoniaMass, mass.first,
                                                             ts.globalMass - MASS_TABLE.ammoniaMass);
                averageAbsError += errors.first;
                averageRelError += errors.second;
                maxAbsError = max(maxAbsError, errors.first);
                maxRelError = max(maxRelError, errors.second);
                minAbsError = min(minAbsError, errors.first);
                minRelError = min(minRelError, errors.second);
                os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                   << "\tsuffix\tammonia loss\tabsolute mass error: " << errors.first <<
                   "\trelative mass error: " << errors.second << "%";
            }
        }
    }
    if (annotatedCnt) {
        averageAbsError /= annotatedCnt;
        averageRelError /= annotatedCnt;
    }
    os << "\n covered " << (double) covered / (ts.sequenceLength - 1) * 100
       << "% of peptide links\ntotal mass: " << ts.globalMass << "\naverage absolute error: " << averageAbsError
       << "\taverage relative error: " << averageRelError << "%\nmax absolute error: " << maxAbsError
       << "\t min absolute error: " << minAbsError <<
       "\nmax relative error: " << maxRelError << "%\t min relative error: " << minRelError << "%\n\n";
}


void printAnnotatedPicks(const tsvParser &p, xmlParser &parser) {
    ofstream os("results.txt");
    for (auto &elem: p.results) {
        string sequence = elem.first;
        os << sequence << "\n";
        TheorySpectra ts = generateMasses(sequence);
        for (auto &spectra: elem.second) {
            showAnnotatedPicks(os, ts, spectra, parser);
        }
    }
    os.close();
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
    parser.parseXML("120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007_msdeconv.msalign");
    tsvParser p;
    p.parseTSV("../120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.tsv");
    printAnnotatedPicks(p, parser);
    return 0;
}