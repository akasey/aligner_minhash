//
// Created by Akash Shrestha on 5/10/18.
//

#include "input_reader.h"


InputReader::InputReader(std::string filename) {
    inputFilename = filename;
    inStream.open(inputFilename, std::ios::in);
    if (!inStream) {
        LOG(ERROR) << "InputReader: Couldn't open " << inputFilename;
        throw InputReaderException("Couldn't open file: " + inputFilename);
    }
}

InputReader::~InputReader(){
    inStream.close();
}