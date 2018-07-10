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
}
