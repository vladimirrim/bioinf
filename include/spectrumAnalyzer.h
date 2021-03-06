#pragma once

#include "massGeneration.h"
#include "xmlParser.h"
#include "tsvParser.h"

class SpectrumAnalyzer {
public:
    static void printAnnotatedPicks(const tsvParser &p, xmlParser &parser, double epsilon, ActivationType filter);

    static void
    printNeutralAndPlusHDifference(const tsvParser &p, xmlParser &parser, const massTable &t, double epsilon);

private:
    static ActivationType filter;
    static const massTable MASS_TABLE;
    static std::vector<std::pair<double, int>> massErrors;
    static std::vector<std::pair<double, int>> massDifferences;
    static double averageRelError, averageAbsError, maxRelError, maxAbsError, minRelError, minAbsError;

    static std::pair<double, double> recalculateError(double theoryMass, double mass);

    static void recalculateMassErrors(double curError);

    static void recalculateMassDifferences(double curDifference, std::vector<std::pair<double, int>> &results);

    static double calculateAnnotatedPicks(const MSResult &spectra, TheorySpectra &ts,
                                          xmlParser &parser, double epsilon, std::ostream &os, bool printData = false);
};