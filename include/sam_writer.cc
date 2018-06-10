//
// Created by Akash Shrestha on 6/9/18.
//

#include "sam_writer.h"

SamWriter::SamWriter(Mode m, std::string filename) : mode(m) {
    file.open(filename);
}

void SamWriter::writeHeaders(std::string wholeCommand, FastaMM &fastamm) {
    std::vector<SamWriter::Header> headers;
    // Add program name to headers
    SamWriter::Header header({"@PG", "ID:"+std::string(PROGRAM_NAME), "CL:"+wholeCommand});
    headers.push_back(header);

    std::map<std::string, int> contigLengthMeta = fastamm.getGenomeLengthMeta();
    for (std::map<std::string, int>::iterator itr=contigLengthMeta.begin(); itr != contigLengthMeta.end(); itr++) {
        std::string serialname = split(itr->first, " ")[0];
        SamWriter::Header header({"@SQ", "SN:"+serialname, "LN:"+std::to_string(itr->second)});
        headers.push_back(header);
    }

    writeHeaders(headers);
}

void SamWriter::writeHeaders(std::vector <SamWriter::Header> &headers) {
    std::for_each(headers.begin(), headers.end(), [this](SamWriter::Header & head){
        head.write(file);
    });
}