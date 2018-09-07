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
#include <mutex>
#include <unordered_map>

#include "common.h"
#include "kmer.h"
#include "murmur3.h"
#include "fixed_array.h"

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
    static const int numHashes = 100;
    static const int LSH_bandSize = 2; // == numHashes/totalBands; choose totalBands such that it's int
    int totalBands = numHashes/LSH_bandSize;
    Kmer MAX_UINT;


#if MINHASH_STORAGE
    std::map<DocID ,Kmer *> minhashStorage;
#endif
    std::vector< std::shared_ptr<std::map<BandhashVar, std::shared_ptr<std::set<DocID > > > > > index;
    std::vector< std::shared_ptr<std::unordered_map<BandhashVar, std::shared_ptr<FixedArray<DocID> > > > > readOnlyIndex;
    static std::mutex mutex;

    void init();
    void deinit();
    std::map<Kmer,int> frequencifyShingles(std::shared_ptr<Kmer> shingles, int len);
    std::shared_ptr<Kmer> computeMinHash(std::map<Kmer,int> shinglesWithFreq);
    std::shared_ptr<BandhashVar> computeBandHash(std::shared_ptr<Kmer> minhash);
    void writeOneIndexToFile(FILE *stream, std::vector<std::shared_ptr<std::map<BandhashVar, std::shared_ptr<std::set<DocID > > > > > &index);
    void loadOneIndexFromFile(FILE *stream, std::vector<std::shared_ptr<std::unordered_map<BandhashVar, std::shared_ptr<FixedArray<DocID > > > > > &r_index, int tb);
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
    void addDocument(DocID id, std::shared_ptr<Kmer> document, int totalShingles);
    std::set<Neighbour> findNeighbours(std::string sequence);
    std::set<Neighbour> findNeighbours(std::shared_ptr<Kmer> doc, int totalShingles);

    void serialize(FILE *stream);
    void deserialize(FILE *stream);
    void deserialize(std::string indexLocation);

    void compareTest(Minhash &second);
};


#endif //ALIGNER_MINHASH_MINHASH_H
