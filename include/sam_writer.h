//
// Created by Akash Shrestha on 6/9/18.
//

#ifndef ALIGNER_MINHASH_SAMWRITER_H
#define ALIGNER_MINHASH_SAMWRITER_H

#include <iostream>
#include "fastamm.h"
#include "ssw_cpp.h"

class SamWriter {
public:
    enum Mode {PAIRED, SINGLE};
    struct Header {
        Header(std::vector<std::string> headers) {
            std::copy(headers.begin(), headers.end(), std::back_inserter(content));
        }
        std::vector<std::string> content;
        void write(std::ofstream &file) {
            for (std::string str : content) {
                file << str << "\t";
            }
            file << std::endl;
        }
    };
    struct Alignment {
        std::string qname;
        uint16_t flag;
        std::string rname;
        uint32_t pos;
        uint8_t mapq;
        std::string cigar;
        std::string rnext;
        uint32_t pnext;
        int32_t tlen;
        std::string seq;
        std::string qual;
        // may be meta
    };

private:
    Mode mode;
    std::ofstream file;
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;

    void writeHeaders(std::vector<Header> &headers);

public:

    SamWriter(Mode m, std::string filename);
    ~SamWriter() {
        file.close();
    }

    void writeHeaders(std::string commandInvoked, FastaMM &fastamm);

};


#endif //ALIGNER_MINHASH_SAMWRITER_H
