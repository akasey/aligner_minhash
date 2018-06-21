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
        segments[segment.segID] = segment;
    }
}

void IndexerJobParser::prepareClassificationJob(FastaMM *fasta) {
    classificationJob.clear();
    for (std::map<int, Segment>::iterator itr=segments.begin(); itr!=segments.end(); itr++) {
        Segment segment = itr->second;
        std::string * sequence = fasta->getGenomePart(segment.key);
        if (!sequence)
            throw ElementNotFoundException(segment.key);
        std::string nucSegment = sequence->substr(segment.start, (segment.end-segment.start) );
        classificationJob[segment.segID] = std::pair<int, std::string>(segment.start, nucSegment);
    }
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

IndexerJobParser::Iterator IndexerJobParser::makeJobIterator() {
    Iterator iterator(&classificationJob);
    return iterator;
};

std::pair<int, std::string> * IndexerJobParser::getSegmentForID(int id) {
    std::map<int, std::pair<int, std::string>>::iterator hit = classificationJob.find(id);
    if (hit == classificationJob.end())
        return NULL;
    return &(hit->second);
};


IndexerJobParser::Iterator::Iterator(std::map<int, std::pair<int, std::string>> *jobs) {
    container = jobs;
    itrPointer = container->begin();
}
bool IndexerJobParser::Iterator::hasNext() {
    return itrPointer!= container->end();
}

std::pair<const int, std::pair<int, std::string>> IndexerJobParser::Iterator::next() {
    std::pair<const int, std::pair<int, std::string>> element = *(itrPointer);
    itrPointer++;
    return element;

};

std::string IndexerJobParser::segIdToKey(int segId) {
    std::map<int, IndexerJobParser::Segment>::iterator found = segments.find(segId);
    if (found != segments.end()) {
        return found->second.key;
    }
    return "";
}