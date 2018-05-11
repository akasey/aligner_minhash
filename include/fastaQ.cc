//
// Created by Akash Shrestha on 5/10/18.
//

#include "fastaQ.h"

bool FastQ::hasNext() {
    int len = inStream.tellg();
    std::string line;
    int linenum = 0;
    for (;linenum<4 && !inStream.eof();linenum++)
        std::getline(inStream, line);
    inStream.seekg(len ,std::ios_base::beg);
    return (linenum==4);
}

InputRead FastQ::next() {
    int linenum = 0;
    std::string line;
    InputRead read;
    while (linenum<4 && std::getline(inStream, line)) {
        switch (linenum) {
            case 0:
                read.key = trim(line);
                break;
            case 1:
                read.sequence = trim(line);
                break;
        }
        linenum++;
    }
    if (linenum !=4 ) {
        LOG(ERROR) << "Fasta file isn't valid";
        throw InputReaderException("Number of lines read: " + std::to_string(linenum) + " is_eof? " + std::to_string(inStream.eof()));
    }
    return read;
}