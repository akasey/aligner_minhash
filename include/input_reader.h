//
// Created by Akash Shrestha on 5/10/18.
//

#ifndef ALIGNER_MINHASH_INPUT_READER_H
#define ALIGNER_MINHASH_INPUT_READER_H

#include "common.h"
#include <fstream>

struct InputRead {
    std::string key;
    std::string sequence;

    InputRead() {

    }

    InputRead(const InputRead &copy) {
        key = copy.key;
        sequence = copy.sequence;
    }
};

struct InputReaderException : public std::exception {
    std::string message;
    InputReaderException(std::string msg) : message(msg) {}
    const char * what () const throw () {
        std::string m = "InputReaderExcpt::" + message;
        return m.c_str();
    }
};

class InputReader {
protected:
    std::string inputFilename;
    std::ifstream inStream;
public:
    InputReader(std::string filename);
    ~InputReader();
    virtual bool hasNext() = 0;
    virtual InputRead next() = 0;
};
#endif //ALIGNER_MINHASH_INPUT_READER_H
