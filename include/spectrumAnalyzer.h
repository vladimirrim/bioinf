#pragma once

#include "massGeneration.h"
#include "xmlParser.h"
#include "tsvParser.h"

class SpectrumAnalyzer{
public:
    static void showAnnotatedPicks(std::ofstream &os, TheorySpectra &ts, const MSResult &spectra, xmlParser &parser, double epsilon);
    static void printAnnotatedPicks(const tsvParser &p, xmlParser &parser, double epsilon,bool filterHCD);
    static void printNeutralAndPlusHDifference(const tsvParser &p, xmlParser &parser, const massTable &t);
private:
    static const massTable MASS_TABLE;
    static double calculateAnnotatedPicks(const MSResult &spectra, const std::vector<MSResult> &elem, TheorySpectra &ts,
                                   xmlParser &parser);
};