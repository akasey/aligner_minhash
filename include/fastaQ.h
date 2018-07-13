//
// Created by Akash Shrestha on 5/10/18.
//

#ifndef ALIGNER_MINHASH_FASTAQ_H
#define ALIGNER_MINHASH_FASTAQ_H

#include "input_reader.h"

/**
 * FastQ is for reading in fastq files.
 */
class FastQ : public InputReader {
public:
    FastQ(std::string filename) : InputReader(filename) {};
    bool hasNext();
    InputRead next();
};

class PairedFastQ {
private:
    FastQ fastq1, fastq2;
public:
    PairedFastQ(std::string files[2]);
    bool hasNext();
    void next(InputRead returnReads[2]);

};

#endif //ALIGNER_MINHASH_FASTAQ_H
