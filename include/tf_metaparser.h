//
// Created by Akash Shrestha on 5/11/18.
//

#ifndef ALIGNER_MINHASH_TF_METAPARSER_H
#define ALIGNER_MINHASH_TF_METAPARSER_H

#include "common.h"
#include <fstream>

class TF_MetaParser {
private:
    std::map<std::string, std::string> meta;
    void init(std::string filename);

public:
    TF_MetaParser(std::string filename);
    std::string& operator[] (std::string key);
    int getInt(std::string key);
};


#endif //ALIGNER_MINHASH_TF_METAPARSER_H
