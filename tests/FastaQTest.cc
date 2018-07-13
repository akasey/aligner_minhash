//
// Created by Akash Shrestha on 5/10/18.
//

#include <iostream>
#include "../include/fastaQ.h"

int main() {
    FastQ reader("../tests/test-inputs/out-subset.bwa.read1.fastq");
//    FastQ reader("/Users/akash/PycharmProjects/aligner/sample_classification_run/out.bwa.read1.fastq");
    while(reader.hasNext()) {
        for (int i=0; i<4 && reader.hasNext(); i++ ) {
            InputRead read = reader.next();
            std::cout << i << " " << read.key << std::endl;
            std::cout << read.sequence << std::endl;
        }
    }

    std::cout << "Testing paired reads" << std::endl;
    std::string inputFiles[2] = { "/Users/akash/PycharmProjects/aligner/sample_classification_run/pairedReads/out.bwa.read1.fastq",
                                  "/Users/akash/PycharmProjects/aligner/sample_classification_run/pairedReads/out.bwa.read2.fastq" };
    PairedFastQ pairedReader(inputFiles);
    int counter = 0;
    while (pairedReader.hasNext() && counter++ < 10) {
        InputRead read[2];
        pairedReader.next(read);
        std::cout << read[0].key << ":" << read[0].sequence << " " << read[1].key << ":" << read[1].sequence << std::endl;
    }
}
