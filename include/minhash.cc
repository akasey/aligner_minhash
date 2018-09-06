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
        index.push_back(std::shared_ptr< std::map<BandhashVar, std::shared_ptr<std::set<DocID> > > >(new std::map<BandhashVar, std::shared_ptr<std::set<DocID> > >()));
    }
}

void Minhash::deinit() {
#if MINHASH_STORAGE
    for (std::map<DocID , Kmer *>::iterator it=minhashStorage.begin(); it!=minhashStorage.end(); it++) {
        delete [] it->second;
    }
#endif
}

std::map<Kmer,int> Minhash::frequencifyShingles(std::shared_ptr<Kmer> shingles, int len) {
    std::map<Kmer, int> count;
    for (int i=0; i<len; i++) {
        Kmer shingle = shingles.get()[i];
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

std::shared_ptr<Kmer> Minhash::computeMinHash(std::map<Kmer,int> shinglesWithFreq) {
    std::shared_ptr<Kmer> hashes( new Kmer[numHashes], [](Kmer* p){ delete [] p; } );
    Kmer bestMinimum[numHashes];
    for (int i=0; i< numHashes; i++)
        bestMinimum[i] = MAX_UINT;

    for (std::map<Kmer,int>::iterator itr = shinglesWithFreq.begin(); itr!=shinglesWithFreq.end(); ++itr) {
        // initial shift value
        Kmer x = itr->first;
        srand(x);
        for (int j=0; j<numHashes; j++){
            for (int k=0; k<itr->second; k++){
                x = rand() % MAX_UINT;
                if (x < bestMinimum[j]) {
                    hashes.get()[j] = itr->first;
                    bestMinimum[j] = x;
                }
            }
        }
    }

    return hashes;
}


std::shared_ptr<BandhashVar> Minhash::computeBandHash(std::shared_ptr<Kmer> minhash) {
    std::shared_ptr<BandhashVar> bandhashes( new BandhashVar[totalBands] , [](BandhashVar* p){ delete[] p; });
    int counter = 0;

    for (int i=0; i<numHashes; i+=LSH_bandSize) {
        Kmer mh[LSH_bandSize];
        for (int j=0; j<LSH_bandSize; j++) {
            mh[j]=minhash.get()[i+j];
        }
        bandhashes.get()[counter] = hashBand(mh,LSH_bandSize);
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
    std::shared_ptr<Kmer> kmers( encodeWindow(sequence, &totalShingles) , [ ](Kmer* p) { delete [ ] p; } );
    addDocument(id, kmers, totalShingles);
}

void Minhash::addDocument(DocID id, std::shared_ptr<Kmer> shingles, int totalShingles) {
    std::map<Kmer,int> shinglesWithFreq = frequencifyShingles(shingles, totalShingles);

    mutex.lock();
    std::shared_ptr<Kmer> minhash = computeMinHash(shinglesWithFreq);
    mutex.unlock();
    std::shared_ptr<BandhashVar> bandhashes = computeBandHash(minhash);

    for(int i=0; i<totalBands; i++) {
        BandhashVar bandHash = bandhashes.get()[i];
        std::map<BandhashVar, std::shared_ptr<std::set<DocID> > >::iterator found = index[i]->find(bandHash);
        if(found == index[i]->end()) {
            std::shared_ptr<std::set<DocID> > a;
            a.get()->insert(id);
            (*index[i])[bandHash] = a;
        }
        else {
            (*index[i])[bandHash]->insert(id);
        }
    }

#if MINHASH_STORAGE
    minhashStorage[id] = minhash;
#endif
}

std::set<Minhash::Neighbour> Minhash::findNeighbours(std::string sequence) {
    int totalShingles = 0;
    std::shared_ptr<Kmer> shingles( encodeWindow(sequence, &totalShingles), [](Kmer* p) { delete [] p; } );
    return findNeighbours(shingles, totalShingles);
}


std::set<Minhash::Neighbour> Minhash::findNeighbours(std::shared_ptr<Kmer> shingles, int totalShingles) {
    std::map<Kmer,int> shinglesWithFreq = frequencifyShingles(shingles, totalShingles);
    std::shared_ptr<Kmer> minhash = computeMinHash(shinglesWithFreq);
    std::shared_ptr<BandhashVar> bandhashes = computeBandHash(minhash);

    std::map<int, int> neighbourhood; // record, hits
    for (int i=0; i<totalBands; i++) {
        BandhashVar bandHash = bandhashes.get()[i];
        std::map<BandhashVar, std::shared_ptr<std::vector<DocID> > >::iterator found = readOnlyIndex[i]->find(bandHash);
        if(found != readOnlyIndex[i]->end()) {
            std::vector<DocID> *hashBand = found->second.get();
//            for (std::set<DocID>::iterator it = hashBand->begin(); it != hashBand->end(); it++) {
            for (const DocID elem : *hashBand) {
                const DocID *it = &elem;
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

    return neighbours;
}


void Minhash::writeOneIndexToFile(FILE *stream, std::vector<std::shared_ptr<std::map<BandhashVar, std::shared_ptr<std::set<DocID > > > > > &index) {
    for (int i=0; i < totalBands; i++) {
        std::map<BandhashVar, std::shared_ptr<std::set<DocID> > > *eachBand = index[i].get();
        int size = eachBand->size();
        write_in_file((void *) &size, sizeof(int), 1, stream);
        for (std::map<BandhashVar, std::shared_ptr<std::set<DocID> > >::iterator it=eachBand->begin(); it!=eachBand->end(); it++) {
            BandhashVar key = it->first;
            std::set<DocID > *value = it->second.get();
            write_in_file((void *) &key, sizeof(BandhashVar), 1, stream);
            int sizeOfVector = value->size();
            write_in_file((void *) &sizeOfVector, sizeof(int), 1, stream);
            for(std::set<DocID>::iterator itr=value->begin(); itr!=value->end(); itr++) {
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

void Minhash::loadOneIndexFromFile(FILE *stream, std::vector<std::shared_ptr<std::map<BandhashVar, std::shared_ptr<std::vector<DocID > > > > > &index, int tb) {
    index.clear();

    for (int i=0; i<tb; i++) {
        std::shared_ptr<std::map<BandhashVar, std::shared_ptr<std::vector<DocID> > > > mapp =  std::shared_ptr<std::map<BandhashVar, std::shared_ptr<std::vector<DocID> > > > (new std::map<BandhashVar, std::shared_ptr<std::vector<DocID> > > ());
        index.push_back( mapp );
        int sizeOfMap = 0;
        read_from_file((void *) &sizeOfMap, sizeof(int), 1, stream);
        for(int j=0; j<sizeOfMap; j++) {
            BandhashVar key;
            read_from_file((void *) &key, sizeof(BandhashVar), 1, stream);
            int sizeOfValue = 0;
            read_from_file((void *) &sizeOfValue, sizeof(int), 1, stream);

            std::shared_ptr<std::vector<DocID> > sett = std::shared_ptr<std::vector<DocID> >(new std::vector<DocID> ());
            (*index[i])[key] = sett;

            for (int k=0; k<sizeOfValue; k++) {
                DocID doc;
                read_from_file((void *) &doc, sizeof(DocID), 1, stream);
                (*index[i])[key].get()->push_back(doc);
            }
        }
    }
}
void Minhash::deserialize(FILE *stream) {
    int tb = 0;
    read_from_file((void *) &tb, sizeof(int), 1, stream);
    totalBands = tb;
    loadOneIndexFromFile(stream, readOnlyIndex, tb);
    LOG(DEBUG) << "Deserialization complete..";
}

void Minhash::deserialize(std::string indexLocation) {
    filename = indexLocation;
    FILE *stream = fopen(indexLocation.c_str(), "rb");
    this->deserialize(stream);
}

void Minhash::compareTest(Minhash &second) {
    for(std::vector<std::shared_ptr<std::map<BandhashVar, std::shared_ptr<std::set<DocID> > > > >::iterator
                itr1 = index.begin(), itr2 = second.index.begin(); itr1 != index.end(); itr1++, itr2++) {
        std::map<BandhashVar, std::shared_ptr<std::set<DocID> > > * first = itr1->get();
        std::map<BandhashVar, std::shared_ptr<std::set<DocID> > > * second = itr2->get();

        for (std::map<BandhashVar, std::shared_ptr<std::set<DocID> > >::iterator itrr = first->begin(); itrr!=first->end(); itrr++) {
            BandhashVar key = itrr->first;
            std::set<DocID > * fromsecond = (*second)[key].get();
            LOG(DEBUG) << (*(itrr->second) == *fromsecond ? "Equal sets" : "Not equal sets");
        }
    }
}
