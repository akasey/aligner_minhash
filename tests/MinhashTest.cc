//
// Created by Akash Shrestha on 5/7/18.
//

#include <iostream>
#include "../include/ThreadPool.h"
#include "../include/minhash.h"

/*
std::map<int, std::string> makeSampleWindows() {
    std::string sequence = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
            "TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
            "TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
            "ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
            "GAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGATGGCAGGTTTCACCG"
            "CCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGC"
            "TGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGT"
            "CAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAGCTTTCCTACTTCGGCG";
    return makeSlidingWindow(sequence, 0, 200, 5);
};


int main(int argc, char *argv[]) {
    Minhash mh;
    std::map<int,std::string> windows = makeSampleWindows();
    for (std::map<int, std::string>::iterator itr=windows.begin(); itr!=windows.end(); itr++) {
        mh.addDocument(itr->first, itr->second);
    }
    FILE *fout = fopen("/Users/akash/ClionProjects/aligner_minhash/tests/test-outputs/minhash.index", "wb");
    mh.serialize(fout);
    fclose(fout);

    Minhash dh;
    FILE *fin = fopen("/Users/akash/ClionProjects/aligner_minhash/tests/test-outputs/minhash.index", "rb");
    dh.deserialize(fin);

    dh.compareTest(mh);
}
*/


void readIndices(std::vector<std::string> mhIndexLocations, int nThreads, std::map<std::string, Minhash *> *mhIndices) {
    ThreadPool threadPool(nThreads);
    for (int i = 0; i < mhIndexLocations.size(); i++) {
        std::string filename = mhIndexLocations[i];
        std::string basefname = basename(filename);
        (*mhIndices)[basefname] = new Minhash();
        threadPool.enqueue([i, filename, &mhIndices, basefname] {
            (*mhIndices)[basefname]->deserialize(filename);
        });
    }
}

int main(int argc, char *argv[]) {
    LOG(INFO) << "Reading in minhash indices...";
//    std::string mhIndexDir = "/Users/akash/PycharmProjects/aligner/sample_classification_run/indices" ;
    std::string mhIndexDir = "/Users/akash/ClionProjects/aligner_minhash/cmake-build-debug/ecoli/indices" ;
    int nThreads = 4;
    std::map<std::string, Minhash *> mhIndices;
    std::vector<std::string> indices = getFilesInDirectory(mhIndexDir, ".mh");
    readIndices(indices, nThreads, &mhIndices);

    while (true) {
        std::string mystr;
        std::cout << "Enter sequence: " << std::endl;
        getline(std::cin, mystr);
        for (int i=0; i<indices.size(); i++) {
            std::string key = "index-" + std::to_string(i) + ".mh";
            std::set<Minhash::Neighbour> neighbours = mhIndices[key]->findNeighbours(mystr);
            if (neighbours.size() > 0) {
                std:: cout << "Neighbours: " << neighbours.size() << " on segment: " << i << std::endl;
            }

        }
        std::cout << "------------------------------" << std::endl;
    }
}