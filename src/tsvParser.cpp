#include <fstream>
#include <map>
#include <cmath>
#include "../include/tsvParser.h"

using namespace std;

void tsvParser::parseTSV(std::string filename) {
    ifstream file(filename);
    string str;
    double threshold = 10e-10;
    getline(file, str);
    while (getline(file, str, '\t')) {
        getline(file, str, '\t');
        MSResult res;
        file >> res.scanId;
        file >> str;
        res.aType = str == "HCD" ? HCD : CID;

        getline(file, str, '\t');
        getline(file, str, '\t');
        getline(file, str, '\t');
        getline(file, str, '\t');
        getline(file, str, '\t');

        file >> res.sequence;
        getline(file, str, '\t');
        getline(file, str, '\t');
        getline(file, str, '\t');
        getline(file, str, '\t');
        getline(file, str, '\t');

        file >> res.eValue;
        if (res.eValue <= threshold) {
            results[res.sequence].push_back(res);
        }

    }
}