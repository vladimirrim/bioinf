#pragma once

#include <string>
#include "xmlParser.h"

struct MSResult{
    int scanId;
    ActivationType aType;
    std::string sequence;
    double eValue;
};

class tsvParser {
public:
    std::map<std::string,std::vector<MSResult>> results;
    void parseTSV(std::string filename);
};
