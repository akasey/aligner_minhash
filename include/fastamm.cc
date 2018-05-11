//
// Created by Akash Shrestha on 5/8/18.
//

#include "fastamm.h"

void FastaMM::readGenome() {
    std::ifstream fin(directory + "/sequence.fasta");
    if (!fin) {
        LOG(ERROR) << "FastaMM.cc: " << directory +"/sequence.fasta not found";
        exit(-1);
    }
    genomeDict.clear();
    std::string line, key;
    while (std::getline(fin, line)) {
        trim(line);
        if (line[0] == '>') {
            key = line;
            genomeDict[key] = "";
        }
        else {
            genomeDict[key] += line;
        }
    }
    fin.close();
}

std::string * FastaMM::getGenomePart(std::string key) {
    std::map<std::string, std::string>::iterator hit = genomeDict.find(key);
    if (hit == genomeDict.end())
        return NULL;
    return &(hit->second);
}