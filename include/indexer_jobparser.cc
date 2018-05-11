//
// Created by Akash Shrestha on 5/8/18.
//

#include "indexer_jobparser.h"


IndexerJobParser::IndexerJobParser(std::string dir) {
    directory = dir;
    std::string filename = directory + "/classify_detail.log";
    std::ifstream fin(filename);
    if (!fin) {
        LOG(ERROR) << "Couldn't open " << filename;
        exit(-1);
    }
    LOG(INFO) << "Reading meta sections from " << filename;
    readMeta(fin);
    LOG(INFO) << "Reading segments from " << filename;
    readSegmentDescription(fin);
    fin.close();
}

void IndexerJobParser::readMeta(std::ifstream &fin) {
    meta.clear();
    std::string line;
    while ( std::getline(fin, line) && line.compare("endMeta")!=0 ) {
        trim(line);
        std::vector<std::string> tokens = split(line, ":");
        if (tokens.size() == 2) {
            meta[trim(tokens[0])] = trim(tokens[1]);
        }
    }
}

void IndexerJobParser::readSegmentDescription(std::ifstream &fin) {
    segments.clear();
    std::string line;
    while(std::getline(fin, line)) {
        trim(line);
        std::vector<std::string> tokens = split(line, "$$$");
        if (tokens.size() != 4) {
            LOG(ERROR) << "Segment description isn't formatted properly. Got: " << line;
            exit(-1);
        }
        Segment segment;
        segment.segID = stoi(tokens[0]);
        segment.key = tokens[1];
        segment.start = stoi(tokens[2]);
        segment.end = stoi(tokens[3]);
        segments.push_back(segment);
    }
}

std::map<int, std::pair<int, std::string>> IndexerJobParser::allClassificationJob(FastaMM *fasta) {
    std::map<int, std::pair<int, std::string>> toReturn;
    for (std::vector<Segment>::iterator itr=segments.begin(); itr!=segments.end(); itr++) {
        std::string * sequence = fasta->getGenomePart(itr->key);
        if (!sequence)
            throw ElementNotFoundException(itr->key);
        std::string segment = sequence->substr(itr->start, (itr->end-itr->start) );
        toReturn[itr->segID] = std::pair<int, std::string>(itr->start, segment);
    }
    return toReturn;
};

int IndexerJobParser::getIntMeta(std::string key) {
    std::map<std::string, std::string>::iterator hit = meta.find(key);
    if (hit == meta.end())
        return NULL;
    return stoi(hit->second);
}

int IndexerJobParser::getWindowLength() {
    return getIntMeta("windowLength");
}

int IndexerJobParser::getNumClasses() {
    return getIntMeta("numClasses");
}

int IndexerJobParser::getK() {
    return getIntMeta("K");
}