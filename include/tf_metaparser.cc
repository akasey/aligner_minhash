//
// Created by Akash Shrestha on 5/11/18.
//

#include "tf_metaparser.h"

void TF_MetaParser::init(std::string filename) {
    std::ifstream inStream(filename);
    if (!inStream) {
        LOG(ERROR) << "TF_MetaParser: Couldn't open file: " << filename;
        exit(-3);
    }

    std::string line;
    while (std::getline(inStream, line)) {
        std::vector<std::string> tokens = split(line, ":");
        if (tokens.size() ==2 ) {
            meta[ trim(tokens[0]) ] = trim(tokens[1]);
        }
    }
}

TF_MetaParser::TF_MetaParser(std::string filename) {
    init(filename);
}

TF_MetaParser::TF_MetaParser(const TF_MetaParser &p2) {
    meta = p2.meta;
}

int TF_MetaParser::getInt(std::string key) {
    std::string hit = (*this)[key];
    return stoi(hit);
}


std::string& TF_MetaParser::operator[] (std::string key) {
    std::map<std::string, std::string>::iterator hit = meta.find(key);
    if (hit == meta.end()) {
        throw TensorflowInferenceException("TF_MetaParser couldn't find key: " + key);
    }
    return hit->second;
}