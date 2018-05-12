#include <iostream>
#include <fstream>
#include "../include/xmlParser.h"
#include "../include/massGeneration.h"
#include "../include/tsvParser.h"
#include "../include/spectrumAnalyzer.h"
#include <cmath>
#include <algorithm>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;

int main() {
    xmlParser parser;
    parser.parseXML("120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007_msdeconv.msalign");
    tsvParser p;
    p.parseTSV("../120706O2c1_LZ-MvD-0297-MabCampth-trypsin_007.tsv");
    double epsilon = 0.3;
    SpectrumAnalyzer::printAnnotatedPicks(p, parser, epsilon, HCD);
    return 0;
}