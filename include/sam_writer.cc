//
// Created by Akash Shrestha on 6/9/18.
//

#include "sam_writer.h"


// Static variables definition
StripedSmithWaterman::Aligner SamWriter::aligner;
StripedSmithWaterman::Filter SamWriter::filter;

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


int SamWriter::alignment(std::string &referenceSegment, std::string &read, SamWriter::Alignment *returnAlignment) {
    int32_t maskLen = strlen(read.c_str())/2;
    maskLen = maskLen < 15 ? 15 : maskLen;
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(read.c_str(), referenceSegment.c_str(), referenceSegment.size(), filter, &alignment, maskLen);

    uint32_t mapq = -4.343 * log(1 - (double)abs(alignment.sw_score - alignment.sw_score_next_best)/(double)alignment.sw_score);
    returnAlignment->mapq = (uint32_t) (mapq + 4.99);
    returnAlignment->cigar = alignment.cigar_string;
    returnAlignment->pos = alignment.ref_begin;

    return (alignment.sw_score/2); // because match score is 2
}