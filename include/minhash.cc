//
// Created by Akash Shrestha on 5/7/18.
//

#include "minhash.h"

#if KMER_K<=8
BandhashVar hashBand(Kmer *minhash, int length) {
    uint32_t out;
    BandhashVar hash = 0;
    for (int i = 0; i < length; i++) {
        Kmer mh = minhash[i];
        MurmurHash3_x86_32(&mh, sizeof(mh), 30, &out);
        hash += out;
    }
    return hash;
}
#elif KMER_K<=16
typedef uint64_t BandhashVar;
BandhashVar hashBand(shingleVar* minhash, int length) {
	BandhashVar hash = 0;
	uint64_t out[2];
	for (int i=0; i<length; i++) {
		Kmer mh = minhash[i];
		MurmurHash3_x64_128 (&mh, sizeof(mh), 30, out);
		hash += out[0];
	}
	return hash;
}
#endif

std::mutex Minhash::mutex;

void Minhash::init() {
    for (int i=0;i<totalBands; i++) {
        index.push_back(new std::map<BandhashVar, std::set<DocID> >());
    }
}

void Minhash::deinit() {
    for (int i=0; i<totalBands; i++) {
        delete index[i];
    }
#if MINHASH_STORAGE
    for (std::map<DocID , Kmer *>::iterator it=minhashStorage.begin(); it!=minhashStorage.end(); it++) {
        delete [] it->second;
    }
#endif
}

std::map<Kmer,int> Minhash::frequencifyShingles(Kmer* shingles, int len) {
    std::map<Kmer, int> count;
    for (int i=0; i<len; i++) {
        Kmer shingle = shingles[i];
        std::map<Kmer,int>::iterator found = count.find(shingle);
        if (found == count.end()){
            count[shingle] = 1;
        }
        else {
            count[shingle] += 1;
        }
    }

    return count;
}

Kmer* Minhash::computeMinHash(std::map<Kmer,int> shinglesWithFreq) {
    Kmer* hashes = new Kmer[numHashes];
    Kmer bestMinimum[numHashes];
    std::fill_n(bestMinimum, numHashes, MAX_UINT);

    for (std::map<Kmer,int>::iterator itr = shinglesWithFreq.begin(); itr!=shinglesWithFreq.end(); ++itr) {
        // initial shift value
        Kmer x = itr->first;
        srand(x);
        for (int j=0; j<numHashes; j++){
            for (int k=0; k<itr->second; k++){
                x = rand() % MAX_UINT;
                if (x < bestMinimum[j]) {
                    hashes[j] = itr->first;
                    bestMinimum[j] = x;
                }
            }
        }
    }

    return hashes;
}


BandhashVar* Minhash::computeBandHash(Kmer* minhash) {
    BandhashVar* bandhashes = new BandhashVar[totalBands];
    int counter = 0;

    for (int i=0; i<numHashes; i+=LSH_bandSize) {
        Kmer mh[LSH_bandSize];
        for (int j=0; j<LSH_bandSize; j++) {
            mh[j]=minhash[i+j];
        }
        bandhashes[counter] = hashBand(mh,LSH_bandSize);
        counter++;
    }
    return bandhashes;
}

float Minhash::jaccard(DocID doc, Kmer* minhash) {
#if MINHASH_STORAGE
    Kmer *doc1S;
    std::map<DocID , Kmer*>::iterator found = minhashStorage.find(doc);
    if (found != minhashStorage.end())
        doc1S = found->second;
    else {
        LOG(ERROR) << doc << " unknown";
        return -1;
    }

    int similarShingles = 0;
    for (int i=0; i<numHashes; i++){
        if(doc1S[i] == minhash[i])
            similarShingles++;
    }

    float similarity = (float)similarShingles/(float)numHashes;
    return similarity;
#else
    return 0.01f;
#endif
}

void Minhash::addDocument(DocID id, std::string sequence) {
    int totalShingles = 0;
    Kmer * kmers = encodeWindow(sequence, &totalShingles);
    Kmer * shingles = new Kmer[totalShingles];
    for (int i=0; i<totalShingles; i++)
        shingles[i] = kmers[i];
    delete [] kmers;
    addDocument(id, shingles, totalShingles);
}

void Minhash::addDocument(DocID id, Kmer *shingles, int totalShingles) {
    std::map<Kmer,int> shinglesWithFreq = frequencifyShingles(shingles, totalShingles);

    mutex.lock();
    Kmer* minhash = computeMinHash(shinglesWithFreq);
    mutex.unlock();
    BandhashVar* bandhashes = computeBandHash(minhash);

    for(int i=0; i<totalBands; i++) {
        std::map<BandhashVar, std::set<DocID > >::iterator found = index[i]->find(bandhashes[i]);
        if(found == index[i]->end()) {
            std::set<DocID > a;
            a.insert(id);
            (*index[i])[bandhashes[i]] = a;
        }
        else {
            (*index[i])[bandhashes[i]].insert(id);
        }
    }

#if MINHASH_STORAGE
    minhashStorage[id] = minhash;
#else
    delete [] minhash;
#endif

    delete [] shingles;
    delete [] bandhashes;
}

std::set<Minhash::Neighbour> Minhash::findNeighbours(std::string sequence) {
    int totalShingles = 0;
    Kmer * kmers = encodeWindow(sequence, &totalShingles);
    Kmer * shingles = new Kmer[totalShingles];
    for (int i=0; i<totalShingles; i++)
        shingles[i] = kmers[i];
    delete [] kmers;
    return findNeighbours(shingles, totalShingles);
}


std::set<Minhash::Neighbour> Minhash::findNeighbours(Kmer* shingles, int totalShingles) {
    std::map<Kmer,int> shinglesWithFreq = frequencifyShingles(shingles, totalShingles);
    Kmer* minhash = computeMinHash(shinglesWithFreq);
    BandhashVar* bandhashes = computeBandHash(minhash);

    std::map<int, int> neighbourhood; // record, hits
    for (int i=0; i<totalBands; i++) {
        std::map<BandhashVar, std::set<DocID > >::iterator found = index[i]->find(bandhashes[i]);
        if(found != index[i]->end()) {
            for (std::set<DocID >::iterator it = found->second.begin();
                 it != found->second.end(); ++it) {
                if (neighbourhood.find(*it) == neighbourhood.end())
                    neighbourhood[*it] = 1;
                else
                    neighbourhood[*it] = neighbourhood[*it] + 1;
            }
        }
    }

    std::set<Minhash::Neighbour> neighbours;
    for (std::map<int, int>::iterator itr = neighbourhood.begin(); itr != neighbourhood.end(); itr++) {
        Minhash::Neighbour rec;
        rec.id = itr->first;
        rec.jaccard = itr->second * 1.0f;
        neighbours.insert(rec);
    }

    delete [] minhash;
    delete [] bandhashes;

    return neighbours;
}


void Minhash::writeOneIndexToFile(FILE *stream, std::vector<std::map<BandhashVar, std::set<DocID > > *> &index) {
    for (int i=0; i < totalBands; i++) {
        std::map<BandhashVar, std::set<DocID > > *eachBand = index[i];
        int size = eachBand->size();
        write_in_file((void *) &size, sizeof(int), 1, stream);
        for (std::map<BandhashVar, std::set<DocID> >::iterator it=eachBand->begin(); it!=eachBand->end(); it++) {
            BandhashVar key = it->first;
            std::set<DocID > value = it->second;
            write_in_file((void *) &key, sizeof(BandhashVar), 1, stream);
            int sizeOfVector = value.size();
            write_in_file((void *) &sizeOfVector, sizeof(int), 1, stream);
            for(std::set<DocID>::iterator itr=value.begin(); itr!=value.end(); itr++) {
                DocID doc = *itr;
                write_in_file((void *) &doc, sizeof(DocID), 1, stream);
            }
        }
    }
}

void Minhash::serialize(FILE *stream) {
    int tb = totalBands;
    write_in_file((void *) &tb, sizeof(int), 1, stream);
    writeOneIndexToFile(stream, index);
    LOG(DEBUG) << "Serialization complete..";
}

void Minhash::loadOneIndexFromFile(FILE *stream, std::vector<std::map<BandhashVar, std::set<DocID > > *> &index, int tb) {
    for (int i=0; i<index.size(); i++) {
        delete index[i];
    }
    index.clear();

    for (int i=0; i<tb; i++) {
        index.push_back(new std::map<BandhashVar, std::set<DocID > >());
        int sizeOfMap = 0;
        read_from_file((void *) &sizeOfMap, sizeof(int), 1, stream);
        for(int j=0; j<sizeOfMap; j++) {
            BandhashVar key;
            read_from_file((void *) &key, sizeof(BandhashVar), 1, stream);
            int sizeOfValue = 0;
            read_from_file((void *) &sizeOfValue, sizeof(int), 1, stream);
            for (int k=0; k<sizeOfValue; k++) {
                DocID doc;
                read_from_file((void *) &doc, sizeof(DocID), 1, stream);
                (*index[i])[key].insert(doc);
            }
        }
    }
}
void Minhash::deserialize(FILE *stream) {
    int tb = 0;
    read_from_file((void *) &tb, sizeof(int), 1, stream);
    totalBands = tb;
    loadOneIndexFromFile(stream, index, tb);
    LOG(DEBUG) << "Deserialization complete..";
}

void Minhash::deserialize(std::string indexLocation) {
    filename = indexLocation;
    FILE *stream = fopen(indexLocation.c_str(), "rb");
    this->deserialize(stream);
}

void Minhash::compareTest(Minhash &second) {
    for(std::vector<std::map<BandhashVar, std::set<DocID > > *>::iterator
                itr1 = index.begin(), itr2 = second.index.begin(); itr1 != index.end(); itr1++, itr2++) {
        std::map<BandhashVar, std::set<DocID > > * first = *itr1;
        std::map<BandhashVar, std::set<DocID > > * second = *itr2;

        for (std::map<BandhashVar, std::set<DocID > >::iterator itrr = first->begin(); itrr!=first->end(); itrr++) {
            BandhashVar key = itrr->first;
            std::set<DocID > fromsecond = (*second)[key];
            LOG(DEBUG) << (itrr->second == fromsecond ? "Equal sets" : "Not equal sets");
        }
    }
}
