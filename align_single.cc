//
// Created by Akash Shrestha on 7/12/18.
//
#include "include/sam_writer.h"
#include "include/cxxopts.h"
#include "include/common.h"
#include "include/minhash.h"
#include "include/tf_metaparser.h"
#include "include/tensorflow_inference.h"
#include "include/ThreadPool.h"
#include "include/fastaQ.h"
#include "include/sam_writer.h"
#include "include/indexer_jobparser.h"
#include "include/priorityqueue_wrapper.h"


struct ReadsWrapper{
    InputRead *read;
    Kmer *kmer;
    Kmer *revKmer;
    int *totalKmers;
    std::string *reverseRead;
    std::set<std::pair<int, bool> > *predictedSegments;

    void clear() {
        if (read) delete read;
        if (kmer) delete [] kmer;
        if (revKmer) delete [] revKmer;
        if (totalKmers) delete totalKmers;
        if (reverseRead) delete reverseRead;
        if (predictedSegments) delete predictedSegments;
    }
};

inline void prediction(TensorflowInference &inferEngine, std::vector<ReadsWrapper > &readsVector, std::vector< std::pair<Kmer *, int> > &pairs, int loadCount) {
    Tensor tensor = inferEngine.makeTensor(pairs);
    std::vector<std::set<std::pair<int, bool> > > predictions = inferEngine.inference(tensor);
    for (int i=0; i<loadCount; i++) {
        readsVector[i].predictedSegments = new std::set<std::pair<int, bool> >();
        std::set<std::pair<int, bool> > inner = predictions[i];
        for (std::set<std::pair<int, bool> >::iterator mnm=predictions[i].begin(); mnm != predictions[i].end(); mnm++) {
            std::pair<int, bool> lacasito = *mnm;
            readsVector[i].predictedSegments->insert(lacasito);
        }
    }
    predictions.clear();
}

inline bool alignMinhashNeighbour(ReadsWrapper *currentRead,
                                  int &predictedSegment, bool &forwardStrand, IndexerJobParser *refBridge, Minhash::Neighbour &neighbour, int windowLength,
                                  SamWriter::Alignment *retAlignment, int *score) {
    std::pair<int, std::string> *referenceSegment = refBridge->getSegmentForID(predictedSegment);
    NULL_CHECK(referenceSegment, "Reference for segment " + std::to_string(predictedSegment) + " is NULL");
    int start = fmax((int)(neighbour.id) - (int)(referenceSegment->first) - 10, 0);
    int length = fmin(referenceSegment->second.length(), start+20+windowLength) - start;
    std::string partOfReference = referenceSegment->second.substr(start, length);
    std::string queryString = forwardStrand ? currentRead->read->sequence : *(currentRead->reverseRead);
    *score = SamWriter::alignment(partOfReference, queryString, retAlignment);

    bool happy;
    retAlignment->qname = currentRead->read->key;
    std::string segmentName = refBridge->segIdToKey(predictedSegment);
//    if (segmentName.empty()) throw
    retAlignment->rname = split(segmentName, " ")[0];
    retAlignment->pos = retAlignment->pos + start + referenceSegment->first;
#if DEBUG_MODE
    LOG(INFO) << "Alignment score: " << *score << " Happy threshold: " << queryString.length()*0.8;
#endif
    if ( queryString.length()*0.8 <= *score) { // consider mapped
        happy = true;
        retAlignment->flag = retAlignment->flag | (forwardStrand ? 0 : REVERSE_MAPPED);
    }
    else {
        happy = false;
        retAlignment->flag = retAlignment->flag | SEGMENT_UNMAPPED;
    }
    return happy;
}

inline bool tryFirstOutofGiven(ReadsWrapper *currentRead, int &predictedSegment, bool &forwardStrand, IndexerJobParser *refBridge,
                               std::set<Minhash::Neighbour> &neighbours, int windowLength, SamWriter::Alignment *retAlignment, int *score) {
    if (neighbours.size() > 0) {
        Minhash::Neighbour first = *(neighbours.begin());
        neighbours.erase(neighbours.begin());
        alignMinhashNeighbour(currentRead, predictedSegment, forwardStrand, refBridge, first, windowLength, retAlignment, score);
        if ( currentRead->read->sequence.length()*0.8 <= *score ) { // atleast 80% matches
            return true;
        }
    }
    return false; // not happy with this alignment
}

inline void makeUnmapped(ReadsWrapper *currentRead, SamWriter::Alignment &retAlignment) {
    retAlignment.qname = currentRead->read->key;
    retAlignment.flag = retAlignment.flag | SEGMENT_UNMAPPED;
}

void align_single(std::string &fastqFile, int &tfBatchSize, TensorflowInference &inferEngine,
                  std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde,
                  SamWriter &samWriter) {
    FastQ fastq(fastqFile);
    while (fastq.hasNext()) {
        std::vector<ReadsWrapper > readsVector(tfBatchSize);
        std::vector< std::pair<Kmer *, int> > pairs(tfBatchSize);
        int loadCount = 0;
        for (int i=0; i<tfBatchSize && fastq.hasNext(); i++) {
            ReadsWrapper readsWrapper;
            readsWrapper.read = new InputRead(fastq.next());
            int totalKmers = 0;
            readsWrapper.kmer = encodeWindow((readsWrapper.read)->sequence, &totalKmers);
            readsWrapper.totalKmers = new int(totalKmers);
            readsVector[i] = readsWrapper;
            std::pair<Kmer *, int> onePair(readsVector[i].kmer, totalKmers);
            pairs[i] = onePair;
            loadCount++;
        }

        prediction(inferEngine, readsVector, pairs, loadCount);

        for (int i=0; i<loadCount; i++) {
            PriorityQueueWrapper queueWrapper;
            ReadsWrapper *currentRead = &(readsVector[i]);
            SamWriter::Alignment bestAlignment;
            int bestScore = -1 * (currentRead->read->sequence.length());
            int score = 0;

            for (std::set<std::pair<int, bool> >::iterator itrr = currentRead->predictedSegments->begin();
                    itrr != currentRead->predictedSegments->end(); itrr++) {
                std::pair<int, bool> pair = *itrr;
                std::string key = "index-" + std::to_string(pair.first) + ".mh";
                if (pair.first < mhIndices.size()) {
                    int partition = pair.first;

                    std::set<Minhash::Neighbour> posNeighboursCurrentPred = mhIndices[key]->findNeighbours(currentRead->kmer, *(currentRead->totalKmers));
#if DEBUG_MODE
                    LOG(INFO) << "+ve Minhash Neigbours: " << posNeighboursCurrentPred.size() << " in segment: " << std::to_string(pair.first);
#endif
                    SamWriter::Alignment alignment;
                    bool forwardStrand = true;
                    bool happy = tryFirstOutofGiven(currentRead, partition, forwardStrand, &referenceGenomeBrigde, posNeighboursCurrentPred,
                                                    referenceGenomeBrigde.getWindowLength(), &alignment, &score);
                    if (!happy) {
                        queueWrapper.addQueue(pair.first, true, &posNeighboursCurrentPred);
                    }
                    else if (score > bestScore) {
                        bestAlignment = alignment;
                        bestScore = score;
                        continue;
                    }


                    // Also try first of reverse strand
                    std::string reverse = reverseComplement(currentRead->read->sequence);
                    currentRead->reverseRead = new std::string(reverse);
                    int totalKmers = 0;
                    currentRead->revKmer = encodeWindow(reverse, &totalKmers);
                    forwardStrand = false;

                    SamWriter::Alignment negAlignment;
                    std::set<Minhash::Neighbour> negNeighboursCurrentPred = mhIndices[key]->findNeighbours(currentRead->revKmer, *(currentRead->totalKmers));
#if DEBUG_MODE
                    LOG(INFO) << "-ve Minhash Neigbours: " << negNeighboursCurrentPred.size() << " in segment: " << std::to_string(pair.first);
#endif
                    happy = tryFirstOutofGiven(currentRead, partition, forwardStrand, &referenceGenomeBrigde, negNeighboursCurrentPred,
                                               referenceGenomeBrigde.getWindowLength(), &negAlignment, &score);
                    if (!happy) {
                        queueWrapper.addQueue(pair.first, false, &negNeighboursCurrentPred);
                    }
                    else {
                        if (score > bestScore) {
                            bestAlignment = negAlignment;
                            bestScore = score;
                        }
                        continue;
                    }

                }
#if DEBUG_MODE
                else {
                    LOG(ERROR) << currentRead->read->key << " predicted as " << std::to_string(pair.first);
                }
#endif
            }

            while(queueWrapper.hasNext()) {
                int partition; bool forwardStrand; Minhash::Neighbour neighbour;
                queueWrapper.pop(&partition, &forwardStrand, &neighbour);
                SamWriter::Alignment retAlignment; int retScore;
                bool happy = alignMinhashNeighbour(currentRead, partition, forwardStrand, &referenceGenomeBrigde, neighbour, referenceGenomeBrigde.getWindowLength(),
                                                   &retAlignment, &retScore);
                if (happy && retScore > bestScore) {
                    bestAlignment = retAlignment;
                    bestScore = score;
                }
            }
            if ( bestAlignment.qname.compare("*")==0 ) { // means no mapping found
                makeUnmapped(currentRead, bestAlignment);
            }
            samWriter.writeAlignment(bestAlignment);
        }

        LOG(INFO) << loadCount << " finished..";
    }
}