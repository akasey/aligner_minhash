//
// Created by Akash Shrestha on 5/8/18.
//

#ifndef ALIGNER_MINHASH_CLASSIFYJOBPARSER_H
#define ALIGNER_MINHASH_CLASSIFYJOBPARSER_H

#include "common.h"
#include "fastamm.h"
#include <map>
#include <fstream>

class IndexerJobParser {
private:
    struct Segment{
        int segID;
        std::string key;
        int start;
        int end;
    };

    std::map<std::string, std::string> meta;
    std::vector<Segment> segments;
    std::string directory;

    void readMeta(std::ifstream &fin);
    int getIntMeta(std::string key);
    void readSegmentDescription(std::ifstream &fin);

public:
    IndexerJobParser(std::string dir);
    ~IndexerJobParser() {}
    std::map<int, std::pair<int, std::string>> allClassificationJob(FastaMM *fasta);
    int getWindowLength();
    int getNumClasses();
    int getK();

};


#endif //ALIGNER_MINHASH_CLASSIFYJOBPARSER_H
