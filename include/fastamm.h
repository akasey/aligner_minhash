//
// Created by Akash Shrestha on 5/8/18.
//

#ifndef ALIGNER_MINHASH_FASTAMM_H
#define ALIGNER_MINHASH_FASTAMM_H

#include "common.h"
#include <fstream>

/**
 * FastaMM takes a directory with sequence.fasta.
 * It reads the file and stores the content as key:string->value:string.
 * key is the fasta entry header for each entries.
 * value is the nucleotide contigs.
 */
class FastaMM {
private:
    std::string directory;
    std::map<std::string, std::string> genomeDict;
    std::map<std::string, int> genomeLengthMeta;
    void readGenome();

public:
    FastaMM(std::string dir) {
        directory = dir;
        readGenome();
    }
    ~FastaMM() {}
    std::string * getGenomePart(std::string key);
    std::map<std::string, int> getGenomeLengthMeta();



};


#endif //ALIGNER_MINHASH_FASTAMM_H
