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
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        case 'N':
        case 'n':
            return rand() % 4;
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
Kmer * encodeWindow(std::string window, int *totalKmers, int K) {
    *totalKmers = window.length() - K + 1;
    Kmer* retHash = new Kmer[*totalKmers];
    Kmer out;
    for (int i=0; i<*totalKmers; i++) {
        std::string kmer = window.substr(i, K);
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

/**
 * Return reverse complement of given sequence
 * @param sequence
 * @return
 */
std::string reverseComplement(std::string sequence) {
    std::transform(sequence.begin(), sequence.end(), sequence.begin(),
            [](char n) {
                switch(n) {
                    case 'A':
                    case 'a':
                        return 'T';
                    case 'T':
                    case 't':
                        return 'A';
                    case 'G':
                    case 'g':
                        return 'C';
                    case 'C':
                    case 'c':
                        return 'G';
                    default:
                        return n;
                }
            });
    std::reverse(sequence.begin(), sequence.end());
    return sequence;
}

