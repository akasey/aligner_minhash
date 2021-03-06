//
// Created by Akash Shrestha on 5/7/18.
//

#ifndef ALIGNER_MINHASH_MINHASH_H
#define ALIGNER_MINHASH_MINHASH_H

#include <map>
#include <vector>
#include <set>
#include <limits>
#include <algorithm>
#include <iostream>

#include "common.h"
#include "kmer.h"
#include "murmur3.h"

#define MINHASH_STORAGE 0

#if KMER_K<=8
typedef uint32_t BandhashVar;
#elif KMER_K<=16
typedef uint64_t BandhashVar;
#endif
BandhashVar hashBand(Kmer * minhash, int length);

class Minhash {
private:
    typedef uint32_t DocID;
    static const int seed = 20; // murmur3 seed
    static const int numHashes = 400;
    static const int LSH_bandSize = 4; // == numHashes/totalBands; choose totalBands such that it's int
    int totalBands = numHashes/LSH_bandSize;
    Kmer MAX_UINT;


#if MINHASH_STORAGE
    std::map<DocID ,Kmer *> minhashStorage;
#endif
    std::vector<std::map<BandhashVar, std::set<DocID > > *> index;

    void init();
    void deinit();
    std::map<Kmer,int> frequencifyShingles(Kmer* shingles, int len);
    Kmer* computeMinHash(std::map<Kmer,int> shinglesWithFreq);
    BandhashVar * computeBandHash(Kmer* minhash);
    void writeOneIndexToFile(FILE *stream, std::vector<std::map<BandhashVar, std::set<DocID > > *> &index);
    void loadOneIndexFromFile(FILE *stream, std::vector<std::map<BandhashVar, std::set<DocID > > *> &index, int tb);
    float jaccard(uint32_t doc, Kmer* minhash);


public:
    std::string filename;
    Minhash() {
        MAX_UINT = std::numeric_limits<Kmer>::max();
        init();
    }

    ~Minhash(){
        deinit();
    }

    struct Neighbour {
        DocID id;
        float jaccard;

        bool operator<(const Neighbour &rhs) const {
            if (jaccard == rhs.jaccard)
                return id > rhs.id;
            return jaccard > rhs.jaccard;
        }

        bool operator==(const Neighbour &rhs) const {
            return id == rhs.id && jaccard==rhs.jaccard;
        }
    };

    void addDocument(DocID id, std::string sequence);
    void addDocument(DocID id, Kmer * document, int totalShingles);
    std::set<Neighbour> findNeighbours(std::string sequence);
    std::set<Neighbour> findNeighbours(Kmer* doc, int totalShingles);

    void serialize(FILE *stream);
    void deserialize(FILE *stream);
    void deserialize(std::string indexLocation);

    void compareTest(Minhash &second);
};


#endif //ALIGNER_MINHASH_MINHASH_H
