//
// Created by Akash Shrestha on 5/7/18.
//

#include "kmer.h"

extern const Kmer POSSIBLE_KMER_SIZE = pow(4,KMER_K);
/**
 *
 * @param a ∈ {A,C,G,T}
 * @return {0,1,2,3} resp.
 */
static uint8_t numbifyNt(char a) {
    switch(a) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
    }
    return -123;
}

/**
 * Converts given kmer to Base4 encoded numeric value.
 * @param kmer
 * @return
 */
Kmer encodeKmer(std::string kmer) {
    int counter = 0;
    Kmer sum = 0;
    for (std::string::iterator it=kmer.begin(); it!=kmer.end(); ++it) {
        sum += numbifyNt(*it) * pow(4,counter);
        counter++;
    }
    return sum;
}

/**
 * Encodes a string window into its contituent numeric kmers.
 * !!! Not a vector of 4**k, just numberic equivalent of kmers
 * @param window: string sequence
 * @return Numeric encoding of kmers in given window and totalKmers = len(window) -K+1
 */
Kmer * encodeWindow(std::string window, int *totalKmers) {
    *totalKmers = window.length() - KMER_K + 1;
    Kmer* retHash = new Kmer[*totalKmers];
    Kmer out;
    for (int i=0; i<*totalKmers; i++) {
        std::string kmer = window.substr(i, KMER_K);
        out = encodeKmer(kmer);
        retHash[i] = out;
    }
    return retHash;
}

/**
 * Computes sliding windows of size `windowLength` from given `segment` jumping `strides` characters in each window.
 * @return
 */
std::map<int, std::string> makeSlidingWindow(std::string &segment, int startLocation, int windowLength, int strides) {
    std::map<int, std::string> windows;
    int end_coord = segment.length() - windowLength;
    for (int start=0; ; start = std::min(start+strides, end_coord) ) {
        std::string sub = segment.substr(start, windowLength);
        windows[startLocation+start] = sub;
        if (start == end_coord)
            break;
    }
    return windows;
}

