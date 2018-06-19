//
// Created by Akash Shrestha on 5/7/18.
//

#ifndef ALIGNER_MINHASH_KMER_H
#define ALIGNER_MINHASH_KMER_H

#include <string>
#include <cmath>
#include <map>

#define KMER_K 7

#if KMER_K<=8
typedef uint16_t Kmer;
#else
typedef uint32_t Kmer;
#endif

extern const Kmer POSSIBLE_KMER_SIZE;

static uint8_t numbifyNt(char a);
Kmer encodeKmer(std::string kmer);
Kmer * encodeWindow(std::string window, int *totalKmers);
std::map<int, std::string> makeSlidingWindow(std::string &segment, int startLocation, int windowLength, int strides);
std::string reverseComplement(std::string sequence);

#endif //ALIGNER_MINHASH_KMER_H
