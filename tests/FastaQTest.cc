//
// Created by Akash Shrestha on 5/10/18.
//

#include <iostream>
#include "../include/fastaQ.h"

int main() {
    FastQ reader("../tests/test-inputs/out-subset.bwa.read1.fastq");
    while(reader.hasNext()) {
        InputRead read = reader.next();
        std::cout << read.key << std::endl;
        std::cout << read.sequence << std::endl;
    }
}
