//
// Created by Akash Shrestha on 5/10/18.
//

#ifndef ALIGNER_MINHASH_FASTAQ_H
#define ALIGNER_MINHASH_FASTAQ_H

#include "input_reader.h"

class FastQ : public InputReader {
public:
    FastQ(std::string filename) : InputReader(filename) {};
    bool hasNext();
    InputRead next();
};


#endif //ALIGNER_MINHASH_FASTAQ_H
