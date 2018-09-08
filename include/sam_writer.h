//
// Created by Akash Shrestha on 6/9/18.
//

#ifndef ALIGNER_MINHASH_SAMWRITER_H
#define ALIGNER_MINHASH_SAMWRITER_H

#include <iostream>
#include <math.h>
#include <mutex>
#include "fastamm.h"
#include "ssw.h"
#include "ssw_cpp.h"

#define PAIRED_READ 1
#define EACH_ALIGNED 2
#define SEGMENT_UNMAPPED 4
#define REVERSE_MAPPED 16



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
        std::string qname = "*";
        uint16_t flag = 0;
        std::string rname = "*";
        uint32_t pos = 0;
        uint8_t mapq = 255;
        std::string cigar = "*";
        std::string rnext = "*";
        uint32_t pnext = 0;
        int32_t tlen = 0;
        std::string seq = "*";
        std::string qual = "*";
        // may be meta

        void print(std::ostream &out) {
            out << qname << " " << pos << " " << mapq << " " << cigar << std::endl;
        }
    };

private:
    Mode mode;
    std::ofstream file;
    std::mutex mutex;
    static StripedSmithWaterman::Aligner aligner;
    static StripedSmithWaterman::Filter filter;

    void writeHeaders(std::vector<Header> &headers);

public:

    SamWriter(Mode m, std::string filename);
    ~SamWriter() {
        file.close();
    }

    void writeHeaders(std::string commandInvoked, FastaMM &fastamm);
    void writeAlignment(Alignment &alignment);
    static int alignment(std::string &referenceSegment, std::string &read, Alignment *returnAlignment);

};


#endif //ALIGNER_MINHASH_SAMWRITER_H
