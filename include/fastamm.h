//
// Created by Akash Shrestha on 5/8/18.
//

#ifndef ALIGNER_MINHASH_FASTAMM_H
#define ALIGNER_MINHASH_FASTAMM_H

#include "common.h"
#include <fstream>

class FastaMM {
private:
    std::string directory;
    std::map<std::string, std::string> genomeDict;
    void readGenome();

public:
    FastaMM(std::string dir) {
        directory = dir;
        readGenome();
    }
    ~FastaMM() {}
    std::string * getGenomePart(std::string key);



};


#endif //ALIGNER_MINHASH_FASTAMM_H
