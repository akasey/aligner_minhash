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
public:
    class Iterator {
    private:
        std::map<int, std::pair<int, std::string>> *container;
        std::map<int, std::pair<int, std::string>>::iterator itrPointer;
    public:
        Iterator(std::map<int, std::pair<int, std::string>> *jobs);
        bool hasNext();
        std::pair<const int, std::pair<int, std::string>> next();
    };
private:
    struct Segment{
        int segID;
        std::string key;
        int start;
        int end;
    };

    std::map<std::string, std::string> meta;
    std::map<int, Segment> segments; // segmentID to Segments
    std::string directory;
    std::map<int, std::pair<int, std::string>> classificationJob; // segmentID -> (start location, sequence)

    void readMeta(std::ifstream &fin);
    int getIntMeta(std::string key);
    void readSegmentDescription(std::ifstream &fin);

public:
    /**
     * Directory containing "classify_detail.log" file in it.
     * This file has some meta informations followed by the neural network classification job.
     * Neural network classifies a read into different segments. The whole nucleotide is read from @FastaMM.
     * The coordinates of the segments are in classify_detail.log file.
     */
    IndexerJobParser(std::string dir);
    ~IndexerJobParser() {}

    /**
     * Takes @FastaMM and using the coordinates of classification job. It creates mapping of segmentID -> (start location, sequence)
     * @param fasta
     * @return
     */
    //std::map<int, std::pair<int, std::string>>
    void prepareClassificationJob(FastaMM *fasta);
    int getWindowLength();
    int getNumClasses();
    int getK();

    Iterator makeJobIterator();
    std::pair<int, std::string> * getSegmentForID(int id);
    std::string segIdToKey(int segId);

};


#endif //ALIGNER_MINHASH_CLASSIFYJOBPARSER_H
