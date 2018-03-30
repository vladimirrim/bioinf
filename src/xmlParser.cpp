#include <fstream>
#include <sstream>
#include "../include/xmlParser.h"

using namespace std;

void xmlParser::parseXML(string filename) {
    ifstream file(filename);
    string str;
    int id;
    while (getline(file, str)) {
        getline(file, str);
        getline(file, str, '=');
        file >> id;
        Spectra &spectra = spectras[id];
        spectra.scanId = id;
        getline(file, str, '=');
        getline(file, str);
        spectra.aType = str == "HCD" ? HCD : CID;
        getline(file, str, '=');
        file >> spectra.precursorMZ;
        getline(file, str, '=');
        file >> spectra.precursorCharge;
        getline(file, str, '=');
        file >> spectra.precursorMass;
        getline(file, str);
        while(getline(file,str)) {
            if(str == "END IONS")
                break;
            string val;
            stringstream is(str);
            getline(is, val, '\t');
            double mass = stod(val);
            getline(is, val, '\t');
            double intensity = stod(val);
            spectra.massesAndIntensities.emplace_back(mass, intensity);
        }
    }
}