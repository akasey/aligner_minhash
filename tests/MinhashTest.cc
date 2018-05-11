//
// Created by Akash Shrestha on 5/7/18.
//

#include "../include/minhash.h"

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
