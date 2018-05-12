#include <fstream>
#include <algorithm>
#include <utility>
#include <opencv/highgui.h>
#include <opencv2/imgproc.hpp>
#include <opencv/cv.hpp>
#include "../include/spectrumAnalyzer.h"

using namespace std;
using namespace cv;

pair<double, double> calculateError(double theoryMass, double mass) {
    double absError = abs(mass - theoryMass);
    double relError = absError / theoryMass * 1e6;
    return {absError, relError};
}


const massTable SpectrumAnalyzer::MASS_TABLE = massTable();
vector<pair<double, int>> SpectrumAnalyzer::massErrors = vector<pair<double, int>>();
vector<pair<double, int>> SpectrumAnalyzer::massDifferences = vector<pair<double, int>>();
ActivationType SpectrumAnalyzer::filter = HCD;
double SpectrumAnalyzer::averageRelError = 0, SpectrumAnalyzer::averageAbsError = 0,
        SpectrumAnalyzer::maxRelError = 0, SpectrumAnalyzer::maxAbsError = 0,
        SpectrumAnalyzer::minRelError = 1e8, SpectrumAnalyzer::minAbsError = 1;

pair<double, double> SpectrumAnalyzer::recalculateError(double theoryMass, double mass) {
    pair<double, double> errors = calculateError(theoryMass, mass);
    averageAbsError += errors.first;
    averageRelError += errors.second;
    maxAbsError = max(maxAbsError, errors.first);
    maxRelError = max(maxRelError, errors.second);
    minAbsError = min(minAbsError, errors.first);
    minRelError = min(minRelError, errors.second);
    return errors;
}


void addToBin(vector<pair<double, int>> &bins, double newValue) {
    for (int i = 0; i < bins.size() - 1; i++) {
        if (newValue > bins[i].first && newValue < bins[i + 1].first) {
            bins[i].second++;
            break;
        }
    }
}

void SpectrumAnalyzer::recalculateMassErrors(double curError) {
    addToBin(massErrors, curError);
}

void SpectrumAnalyzer::recalculateMassDifferences(double curDifference,
                                                  vector<pair<double, int>> &results) {
    addToBin(massDifferences, curDifference);

    addToBin(results, curDifference);
}

void drawHCDHistogram(string filename, const vector<pair<double, int>> &massDifferences) {
    int hist_w = 1400;
    int hist_h = 670;
    int bin_w = cvRound((double) (hist_w - 30) / massDifferences.size());
    Mat histImage(hist_h, hist_w, CV_8UC3, Scalar(0, 0, 0));
    int XScale = -40;
    for (int i = 0; i < massDifferences.size(); i++) {
        rectangle(histImage, Point(bin_w * i + 35, hist_h - 30),
                  Point(bin_w * i + bin_w + 25, hist_h - 30 - massDifferences[i].second),
                  Scalar(130, 221, 238), -1);
        putText(histImage, (string) (XScale < 0 ? "-" : "") + "0.0" +
                           (abs(XScale) % 10 == 0 ? to_string(abs(XScale) / 10) : "0" + to_string(abs(XScale))),
                Point(bin_w * i + (XScale < 0 ? 15 : 20), hist_h - 10), FONT_HERSHEY_PLAIN,
                0.8, Scalar(255, 255, 255));
        XScale += 4;
    }
    int YScale = 0;
    for (int i = 0; i < 21; i++) {
        putText(histImage, to_string(YScale), Point(10, hist_h - YScale - 30), FONT_HERSHEY_PLAIN, 0.8,
                Scalar(255, 255, 255));
        YScale += (hist_h - 30) / 20;
    }
    imwrite(filename, histImage);
}

void drawCIDHistogram(string filename, const vector<pair<double, int>> &massDifferences) {
    int hist_w = 1200;
    int hist_h = 740;
    int bin_w = cvRound((double) (hist_w - 30) / massDifferences.size());
    Mat histImage(hist_h, hist_w, CV_8UC3, Scalar(0, 0, 0));
    int XScale = -60;
    for (int i = 0; i < massDifferences.size(); i++) {
        rectangle(histImage, Point(bin_w * i + 35, hist_h - 30),
                  Point(bin_w * i + bin_w + 25, hist_h - 30 - 15 * massDifferences[i].second),
                  Scalar(130, 221, 238), -1);
        putText(histImage, (string) (XScale < 0 ? "-" : "") + "0." + (abs(XScale) == 6 ? "0" : "") +
                           (abs(XScale) % 10 == 0 ? to_string(abs(XScale) / 10) : to_string(abs(XScale))),
                Point(bin_w * i + (XScale < 0 ? 15 : 20), hist_h - 10), FONT_HERSHEY_PLAIN,
                0.8, Scalar(255, 255, 255));
        XScale += 6;
    }
    int YScale = 0;
    for (int i = 0; i < 24; i++) {
        putText(histImage, to_string(YScale), Point(10, hist_h - 15 * YScale - 30), FONT_HERSHEY_PLAIN, 0.8,
                Scalar(255, 255, 255));
        YScale += (hist_h - 30) / 345;
    }
    imwrite(filename, histImage);
}

void
SpectrumAnalyzer::printAnnotatedPicks(const tsvParser &p, xmlParser &parser, double epsilon, ActivationType filter) {
    ofstream os("results.txt");
    SpectrumAnalyzer::filter = filter;
    if (filter == HCD) {
        double start = -0.04;
        for (int i = 0; i <= 20; i++) {
            massErrors.emplace_back(start, 0);
            massDifferences.emplace_back(start, 0);
            start += 0.004;
        }
    } else {
        double start = -0.6;
        for (int i = 0; i <= 20; i++) {
            massErrors.emplace_back(start, 0);
            massDifferences.emplace_back(start, 0);
            start += 0.06;
        }
    }
    for (auto &elem: p.results) {
        string sequence = elem.first;
        os << sequence << "\n";
        TheorySpectra ts = generateMasses(sequence);
        for (auto &spectra: elem.second) {
            if (spectra.aType == filter)
                calculateAnnotatedPicks(spectra, ts, parser, epsilon, os, true);
        }
    }
    sort(massErrors.begin(), massErrors.end());
    os << "Exact mass errors:\n";
    for (auto &error:massErrors) {
        os << "Error: " << error.first << "\t Times occurred: " << error.second << "\n";
    }

    sort(massDifferences.begin(), massDifferences.end());
    os << "\nMass differences:\n";
    for (auto &difference:massDifferences) {
        os << "Error: " << difference.first << "\t Times occurred: " << difference.second << "\n";
    }
    os.close();
    if (filter == HCD)
        drawHCDHistogram("HCD.png", massDifferences);
    else
        drawCIDHistogram("CID.png", massDifferences);
}

void SpectrumAnalyzer::printNeutralAndPlusHDifference(const tsvParser &p, xmlParser &parser, const massTable &t,
                                                      double epsilon) {
    ofstream os("difference.txt");
    vector<pair<double, MSResult>> differences;
    for (auto &elem: p.results) {
        string sequence = elem.first;
        TheorySpectra tsNeutral = generateMasses(sequence);
        TheorySpectra tsPlusH = generateMasses(sequence, t.HMass);
        for (auto &spectra: elem.second) {
            double coveredNeutral = calculateAnnotatedPicks(spectra, tsNeutral, parser, epsilon, os);
            double coveredPlusH = calculateAnnotatedPicks(spectra, tsPlusH, parser, epsilon, os);
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
        calculateAnnotatedPicks(spectra, tsNeutral, parser, epsilon, os, true);
        os << "Plus H: \n";
        calculateAnnotatedPicks(spectra, tsPlusH, parser, epsilon, os, true);
    }
    os.close();
}


double
SpectrumAnalyzer::calculateAnnotatedPicks(const MSResult &spectra, TheorySpectra &ts,
                                          xmlParser &parser, double epsilon, ostream &os, bool printData) {
    if (printData)
        os << "\n" << spectra.scanId << "\t" << spectra.eValue << "\t" << (spectra.aType == HCD ? "HCD" : "CID");

    double massEpsilon = 0.001;
    int covered = 0;
    map<string, bool> isCovered;
    map<string, double> annotatedPrefixes;
    auto comp = [](const string &a, const string &b) { return a.length() < b.length(); };
    map<string, double, decltype(comp)> annotatedSuffixes(comp);
    vector<pair<double, int>> prefDifferences;
    vector<pair<double, int>> sufDifferences;
    if (filter == HCD) {
        double start = -0.1;
        for (int i = 0; i <= 20; i++) {
            prefDifferences.emplace_back(start, 0);
            sufDifferences.emplace_back(start, 0);
            start += 0.01;
        }
    } else {
        double start = -1;
        for (int i = 0; i <= 20; i++) {
            prefDifferences.emplace_back(start, 0);
            sufDifferences.emplace_back(start, 0);
            start += 0.1;
        }
    }
    int annotatedCnt = 0;
    averageRelError = 0, averageAbsError = 0, maxRelError = 0, maxAbsError = 0, minRelError = 1e8, minAbsError = 1;
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
                auto curError = (double) (mass.first - prefix.second);
                recalculateMassErrors(curError);
                annotatedPrefixes[prefix.first] = mass.first;
                pair<double, double> errors = recalculateError(prefix.second, mass.first);
                if (printData) {
                    os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                       << "\tprefix\t\tabsolute mass error: "
                       << errors.first <<
                       "\trelative mass error: " << errors.second << "ppm\tExact mass error: "
                       << curError;
                }
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.waterMass) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                auto curError = (double) (mass.first - prefix.second + MASS_TABLE.waterMass);
                recalculateMassErrors(curError);
                annotatedPrefixes[prefix.first] = mass.first + MASS_TABLE.waterMass;
                pair<double, double> errors = recalculateError(prefix.second - MASS_TABLE.waterMass, mass.first);
                if (printData) {
                    os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                       << "\tprefix\twater loss\tabsolute mass error: " << errors.first <<
                       "\trelative mass error: " << errors.second << "ppm\tExact mass error: "
                       << (double) (mass.first - prefix.second + MASS_TABLE.waterMass);
                }
            }

            if (abs(mass.first - prefix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[prefix.first]) {
                    covered++;
                    isCovered[prefix.first] = true;
                    isCovered[ts.linkedPairs[prefix.first]] = true;
                }
                auto curError = (double) (mass.first - prefix.second + MASS_TABLE.ammoniaMass);
                recalculateMassErrors(curError);
                annotatedPrefixes[prefix.first] = mass.first + MASS_TABLE.ammoniaMass;
                pair<double, double> errors = recalculateError(prefix.second - MASS_TABLE.ammoniaMass, mass.first);
                if (printData) {
                    os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                       << "\tprefix\tammonia loss\tabsolute mass error: " << errors.first <<
                       "\trelative mass error: " << errors.second << "ppm\tExact mass error: "
                       << curError;
                }
            }

        }
    }
    for (auto &suffix:ts.suffixes) {
        string name = suffix.first;
        for (auto &mass: parser.spectras[spectra.scanId].massesAndIntensities) {
            if (abs(mass.first - suffix.second) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[suffix.first]) {
                    covered++;
                    isCovered[suffix.first] = true;
                    isCovered[ts.linkedPairs[suffix.first]] = true;
                }
                auto curError = (double) (mass.first - suffix.second);
                recalculateMassErrors(curError);
                annotatedSuffixes[suffix.first] = mass.first;
                pair<double, double> errors = recalculateError(suffix.second, mass.first);
                if (printData) {
                    os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                       << "\tsuffix\t\tabsolute mass error: "
                       << errors.first <<
                       "\trelative mass error: " << errors.second << "ppm\tExact mass error: "
                       << (double) (mass.first - suffix.second);
                }
            }

            if (abs(mass.first - suffix.second + MASS_TABLE.waterMass) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[suffix.first]) {
                    covered++;
                    isCovered[suffix.first] = true;
                    isCovered[ts.linkedPairs[suffix.first]] = true;
                }
                auto curError = (double) (mass.first - suffix.second + MASS_TABLE.waterMass);
                recalculateMassErrors(curError);
                annotatedSuffixes[suffix.first] = mass.first + MASS_TABLE.waterMass;
                pair<double, double> errors = recalculateError(suffix.second - MASS_TABLE.waterMass, mass.first);
                if (printData) {
                    os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                       << "\tsuffix\twater loss\tabsolute mass error: " << errors.first <<
                       "\trelative mass error: " << errors.second << "ppm\tExact mass error: "
                       << curError;
                }
            }

            if (abs(mass.first - suffix.second + MASS_TABLE.ammoniaMass) <= epsilon) {
                annotatedCnt++;
                if (!isCovered[suffix.first]) {
                    covered++;
                    isCovered[suffix.first] = true;
                    isCovered[ts.linkedPairs[suffix.first]] = true;
                }
                auto curError = (double) (mass.first - suffix.second + MASS_TABLE.ammoniaMass);
                recalculateMassErrors(curError);
                annotatedSuffixes[suffix.first] = mass.first + MASS_TABLE.ammoniaMass;
                pair<double, double> errors = recalculateError(suffix.second - MASS_TABLE.ammoniaMass, mass.first);
                if (printData) {
                    os << "\n" << mass.first << "\t" << mass.second << "\t" << name
                       << "\tsuffix\tammonia loss\tabsolute mass error: " << errors.first <<
                       "\trelative mass error: " << errors.second << "ppm\tExact mass error: "
                       << (double) (mass.first - suffix.second + MASS_TABLE.ammoniaMass);
                }
            }
        }
    }
    if (annotatedCnt) {
        averageAbsError /= annotatedCnt;
        averageRelError /= annotatedCnt;
    }
    pair<string, double> prevPick = {"1", 0};
    for (auto &prefix: annotatedPrefixes) {
        if (prevPick.first != "1" && prefix.first.length() - prevPick.first.length() <= 3) {
            double theoryDifference = calculatePeptideMass(prefix.first) - calculatePeptideMass(prevPick.first);
            double actualDifference = prefix.second - prevPick.second;
            recalculateMassDifferences(theoryDifference - actualDifference, prefDifferences);
        }
        prevPick = prefix;
    }
    prevPick = {"1", 0};
    for (auto &suffix: annotatedSuffixes) {
        if (prevPick.first != "1" && suffix.first.length() - prevPick.first.length() <= 3) {
            double theoryDifference = calculatePeptideMass(suffix.first) - calculatePeptideMass(prevPick.first);
            double actualDifference = suffix.second - prevPick.second;
            recalculateMassDifferences(theoryDifference - actualDifference, sufDifferences);
        }
        prevPick = suffix;
    }

    if (printData) {
        os << "\n covered " << (double) covered / (ts.sequenceLength - 1) * 100
           << "% of peptide links\ntotal mass: " << ts.globalMass << "\naverage absolute error: " << averageAbsError
           << "\taverage relative error: " << averageRelError << "ppm\nmax absolute error: " << maxAbsError
           << "\t min absolute error: " << minAbsError <<
           "\nmax relative error: " << maxRelError << "ppm\t min relative error: " << minRelError << "ppm\n";

        os << "\nMass differences:\n";
        os << "Prefixes:\n";
        for (auto &difference:prefDifferences) {
            os << "Difference: " << difference.first << "\t Times occurred: " << difference.second << "\n";
        }
        os << "Suffixes:\n";
        for (auto &difference:sufDifferences) {
            os << "Difference: " << difference.first << "\t Times occurred: " << difference.second << "\n";
        }
        os << "\n";
    }
    return (double) covered / (ts.sequenceLength - 1) * 100;
}

