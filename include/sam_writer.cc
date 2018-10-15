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


int SamWriter::alignment(std::string &referenceSegment, std::string &read, SamWriter::Alignment *returnAlignment, int *num_mismatches) {
    int32_t maskLen = strlen(read.c_str())/2;
    maskLen = maskLen < 15 ? 15 : maskLen;
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(read.c_str(), referenceSegment.c_str(), referenceSegment.size(), filter, &alignment, maskLen);

    uint32_t mapq = -4.343 * log(1 - (double)abs(alignment.sw_score - alignment.sw_score_next_best)/(double)alignment.sw_score);
    returnAlignment->mapq = (uint32_t) (mapq + 4.99);
    returnAlignment->cigar = alignment.cigar_string;
    returnAlignment->pos = alignment.ref_begin + 1;

    returnAlignment->seq = std::string(read);
    *num_mismatches = alignment.mismatches;
/*
    // All below for alignment SEQ....
    char alignSEQ[read.length()+1];
    int alignPtr = 0; // where to fill up in alignSEQ array
    int start = 0; // which index to copy from read. Might have to skip if there is deletion from reference
    for (int i=0; i<alignment.cigar.size(); i++) {
        uint32_t bam = alignment.cigar[i];
        int length = cigar_int_to_len(bam);
        char op = cigar_int_to_op(bam);
        switch (op) {
            case 'M':
            case 'I':
            case 'S':
            case '=':
            case 'X':
                for (int j=0;j<length; j++) {
                    alignSEQ[alignPtr++] = read[start++];
                }
                break;
        }
    }
    alignSEQ[alignPtr] = '\0';
    returnAlignment->seq = std::string(alignSEQ);
    std::cout << returnAlignment->seq.length() << " " << alignPtr << std:: endl;
    // end of alignment SEQ
*/

    return (alignment.sw_score/2); // because match score is 2
}

int SamWriter::alignmentEdlib(std::string &referenceSegment, std::string &read, SamWriter::Alignment *returnAlignment, int *num_mismatches) {
    *num_mismatches = read.length();
    EdlibAlignTask task = EDLIB_TASK_PATH;
    EdlibAlignMode mode = EDLIB_MODE_HW;

    char lowercase[] = {'a','c','g','t'};
    char uppercase[] = {'A', 'C', 'G', 'T'};
    EdlibEqualityPair pairs[4];
    for (int i=0; i<4; i++) {
        EdlibEqualityPair pair;
        pair.first = lowercase[i];
        pair.second = uppercase[i];
        pairs[i] = pair;
    }
    EdlibAlignResult result = edlibAlign(read.c_str(), read.length(), referenceSegment.c_str(), referenceSegment.length(),
                                         edlibNewAlignConfig(-1, mode, task, pairs, 4));
    if (result.status == EDLIB_STATUS_OK) {
        char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
        returnAlignment->mapq = 255; // not available
        returnAlignment-> cigar = std::string(cigar);
        returnAlignment->pos = *(result.startLocations);
        returnAlignment->seq = std::string(read);
        *num_mismatches = result.editDistance;
    }
    edlibFreeAlignResult(result);
    return read.length() - *num_mismatches;
}

void SamWriter::writeAlignment(Alignment &alignment) {
    std::unique_lock<std::mutex> lock(mutex);
    file << alignment.qname << "\t"
        << unsigned(alignment.flag) << "\t"
        << alignment.rname << "\t"
        << unsigned(alignment.pos) << "\t"
        << unsigned(alignment.mapq) << "\t"
        << alignment.cigar << "\t"
        << alignment.rnext << "\t"
        << unsigned(alignment.pnext) << "\t"
        << alignment.tlen << "\t"
        << alignment.seq << "\t"
        << alignment.qual << "\n";
    file.flush();
}