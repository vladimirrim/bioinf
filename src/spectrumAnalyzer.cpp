#include <fstream>
#include <algorithm>
#include "../include/spectrumAnalyzer.h"

using namespace std;


pair<double, double> calculateError(double theoryMass, double mass) {
    double absError = abs(mass - theoryMass);
    double relError = absError / theoryMass * 100;
    return {absError, relError};
}

const massTable SpectrumAnalyzer::MASS_TABLE = massTable();

void
SpectrumAnalyzer::showAnnotatedPicks(std::ofstream &os, TheorySpectra &ts, const MSResult &spectra, xmlParser &parser,
                                     double epsilon) {
    os << "\n" << spectra.scanId << "\t" << spectra.eValue << "\t" << (spectra.aType == HCD ? "HCD" : "CID");
    int covered = 0;
    map<string, bool> isCovered;
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
                pair<double, double> errors = calculateError(prefix.second, mass.first);
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
                pair<double, double> errors = calculateError(prefix.second - MASS_TABLE.waterMass, mass.first);
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
                pair<double, double> errors = calculateError(prefix.second - MASS_TABLE.ammoniaMass, mass.first);
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
                pair<double, double> errors = calculateError(prefix.second, mass.first);
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
                pair<double, double> errors = calculateError(prefix.second - MASS_TABLE.waterMass, mass.first);
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
                pair<double, double> errors = calculateError(prefix.second - MASS_TABLE.ammoniaMass, mass.first);
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

void SpectrumAnalyzer::printAnnotatedPicks(const tsvParser &p, xmlParser &parser, double epsilon, bool filterHCD) {
    ofstream os("results.txt");
    for (auto &elem: p.results) {
        string sequence = elem.first;
        os << sequence << "\n";
        TheorySpectra ts = generateMasses(sequence);
        for (auto &spectra: elem.second) {
            if (spectra.aType == CID || !filterHCD)
                showAnnotatedPicks(os, ts, spectra, parser, epsilon);
        }
    }
    os.close();
}

void SpectrumAnalyzer::printNeutralAndPlusHDifference(const tsvParser &p, xmlParser &parser, const massTable &t) {
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
    double epsilon = 0.02;
    for (int i = 1; i <= 3; i++) {
        MSResult spectra = differences[differences.size() - i].second;
        string sequence = spectra.sequence;
        TheorySpectra tsNeutral = generateMasses(sequence);
        TheorySpectra tsPlusH = generateMasses(sequence, t.HMass);
        os << sequence << "\n";
        os << "difference: " << differences[differences.size() - i].first << "\n";
        os << "Neutral:\n";
        showAnnotatedPicks(os, tsNeutral, spectra, parser, epsilon);
        os << "Plus H: \n";
        showAnnotatedPicks(os, tsPlusH, spectra, parser, epsilon);
    }
    os.close();
}


double
SpectrumAnalyzer::calculateAnnotatedPicks(const MSResult &spectra, const std::vector<MSResult> &elem, TheorySpectra &ts,
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