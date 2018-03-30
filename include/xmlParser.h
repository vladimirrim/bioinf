#pragma once

#include <utility>
#include <vector>
#include <string>
#include <map>

enum ActivationType {
    HCD, CID
};

struct Spectra {
    int scanId;
    ActivationType aType;
    double precursorMass, precursorCharge, precursorMZ;
    std::vector<std::pair<double, double>> massesAndIntensities;
};

class xmlParser {
public:
    std::map<int,Spectra> spectras;
    void parseXML(std::string filename);
};