//
// Created by Akash Shrestha on 7/12/18.
//

#ifndef ALIGNER_ALIGN_SINGLE_H
#define ALIGNER_ALIGN_SINGLE_H

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

#define NBEST 3
#define EXTRA_EXTENSION_SHORT 100

struct ReadsWrapper{
    std::shared_ptr<InputRead> read;
    std::shared_ptr<Kmer> kmer;
    std::shared_ptr<Kmer> revKmer;
    std::shared_ptr<int> totalKmers;
    std::shared_ptr<std::string> reverseRead;
    std::shared_ptr<std::set<std::pair<int, bool> > > predictedSegments;

    ReadsWrapper() {
        std::nullptr_t null;
        read = std::shared_ptr<InputRead>(null);
        kmer = std::shared_ptr<Kmer>(null, [](Kmer* p) { delete [] p; });
        revKmer = std::shared_ptr<Kmer>(null, [](Kmer* p) { delete [] p; });
        totalKmers = std::shared_ptr<int>(null);
        reverseRead = std::shared_ptr<std::string>(null);
        predictedSegments = std::shared_ptr<std::set<std::pair<int, bool> > >(null);
    }
};

inline void prediction(TensorflowInference &inferEngine, std::vector<ReadsWrapper > &readsVector, std::vector< std::pair<std::shared_ptr<Kmer>, int> > &pairs, int loadCount) {
    Tensor tensor = inferEngine.makeTensor(pairs);
    std::vector<std::set<std::pair<int, bool> > > predictions = inferEngine.inference(tensor);
    for (int i=0; i<loadCount; i++) {
        readsVector[i].predictedSegments.reset( new std::set<std::pair<int, bool> >() );
        std::set<std::pair<int, bool> > inner = predictions[i];
        for (std::set<std::pair<int, bool> >::iterator mnm=predictions[i].begin(); mnm != predictions[i].end(); mnm++) {
            std::pair<int, bool> lacasito = *mnm;
            readsVector[i].predictedSegments->insert(lacasito);
        }
    }
    predictions.clear();
}

inline void insertNewAlignment(std::vector<std::pair<SamWriter::Alignment,int>> &bestAlignments, std::pair<SamWriter::Alignment,int> newAlignment, int threshold) {
    bool inserted = false;
    for (int i=0; i<bestAlignments.size(); i++) {
        std::pair<SamWriter::Alignment, int> *each = &bestAlignments[i];
        if (absoluteValue((int)(each->first.pos - newAlignment.first.pos)) <= threshold) {
            if (each->second <= newAlignment.second) {
                bestAlignments[i] = newAlignment;
                inserted = true;
            }
        }
    }
    if (!inserted)
        bestAlignments.push_back(newAlignment);
}

inline bool isNeigbourPreviouslySeen(std::unordered_set<int> &seenNeigbbours, Minhash::Neighbour &neighbour, int distance) {
    // If this neighbour is within half of the distance, then return true
    for (const int &seenNeighbour : seenNeigbbours) {
        if ( absoluteValue(((int)neighbour.id - seenNeighbour)) <= EXTRA_EXTENSION_SHORT/2 ) {
            return true;
        }
    }
    return false;
}

inline bool alignMinhashNeighbour(ReadsWrapper *currentRead,
                                  int &predictedSegment, bool &forwardStrand, IndexerJobParser *refBridge, Minhash::Neighbour &neighbour,
                                  SamWriter::Alignment *retAlignment, int *score) {
    std::pair<int, std::string> *referenceSegment = refBridge->getSegmentForID(predictedSegment);
    NULL_CHECK(referenceSegment, "Reference for segment " + std::to_string(predictedSegment) + " is NULL");
    int windowLength = currentRead->read->sequence.length();
    int start = fmax((int)(neighbour.id) - (int)(referenceSegment->first) - EXTRA_EXTENSION_SHORT, 0);
    int length = fmin(referenceSegment->second.length(), start+EXTRA_EXTENSION_SHORT*2+windowLength) - start;
    int numMismatches = 0;
    std::string partOfReference = referenceSegment->second.substr(start, length);
    std::string queryString = forwardStrand ? currentRead->read->sequence : *(currentRead->reverseRead);
    *score = SamWriter::alignment(partOfReference, queryString, retAlignment, &numMismatches);

    bool happy;
    retAlignment->qname = currentRead->read->key;
    std::string segmentName = refBridge->segIdToKey(predictedSegment);
//    if (segmentName.empty()) throw
    retAlignment->rname = split(segmentName, " ")[0];
    retAlignment->pos = retAlignment->pos + start + referenceSegment->first;
#if DEBUG_MODE
    LOG(INFO) << "Alignment score: " << *score << " Happy threshold: " << queryString.length()*0.8;
#endif
    if ( queryString.length()*0.8 <= *score) { // atleast 80% matches // consider mapped
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
                               std::set<Minhash::Neighbour> &neighbours, SamWriter::Alignment *retAlignment, int *score, std::unordered_set<int> &seenNeighbours) {
    if (neighbours.size() > 0) {
        Minhash::Neighbour first;
        do {
            first = *(neighbours.begin());
            neighbours.erase(neighbours.begin());
        } while (neighbours.size()>0 && isNeigbourPreviouslySeen(seenNeighbours, first, refBridge->getWindowLength()));
        // atleast 80% matches
        return alignMinhashNeighbour(currentRead, predictedSegment, forwardStrand, refBridge, first, retAlignment, score);
    }
    return false; // not happy with this alignment
}

inline void makeUnmapped(ReadsWrapper *currentRead, SamWriter::Alignment &retAlignment) {
    retAlignment.qname = currentRead->read->key;
    retAlignment.flag = retAlignment.flag | SEGMENT_UNMAPPED;
}


/****
 *
 * @NOTE to myself.
 * Dude you're indexing both positive window and negative windows into minhash. - Oh no, we're indexing positive windows only
 * You are predicting with neural network just once... You could do the same with minhash search as well.
 * The things that you put into priority queue are just redundant.
 *
 * Now if you search minhash just once, you won't know which strand you're dealing with.
 * So you'll end up trying both strand. Is it expensive this way?
*****/
std::vector<std::pair<SamWriter::Alignment, int>> alignOne(ReadsWrapper eachRead, std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde, int nBest) {
    ReadsWrapper *currentRead = &eachRead;
    PriorityQueueWrapper queueWrapper;
    std::vector<std::pair<SamWriter::Alignment,int>> bestAlignments;
    std::unordered_set<int> seenNeighbours; // list of neighbourhood where bestAlignments were already found. So we won't repeat in same territory
    int bestScore = -1 * (currentRead->read->sequence.length());
    float bestNeighbourJaccardScore = 0.0; // best jaccard score of the minhash neighbour yet.. If the next neighbour we are searching is so worse we terminate search
    int score = 0;

#if DEBUG_MODE
    LOG(INFO) << "Predicted size :" << currentRead->predictedSegments->size();
#endif

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
            bestNeighbourJaccardScore = std::max(bestNeighbourJaccardScore,posNeighboursCurrentPred.begin()->jaccard);
            SamWriter::Alignment alignment;
            bool forwardStrand = true;
            std::pair<int, bool> firstResult = { posNeighboursCurrentPred.begin()->id, forwardStrand};
            bool happy = tryFirstOutofGiven(currentRead, partition, forwardStrand, &referenceGenomeBrigde, posNeighboursCurrentPred,
                                            &alignment, &score, seenNeighbours);

            if (happy) {
                insertNewAlignment(bestAlignments, {alignment, score}, EXTRA_EXTENSION_SHORT);
                seenNeighbours.insert(firstResult.first);
                if (score > bestScore)                 bestScore = score;
                if (bestAlignments.size() >= nBest)    continue;
            }
            queueWrapper.addQueue(pair.first, true, &posNeighboursCurrentPred);


            // Also try first of reverse strand
            std::string reverse = reverseComplement(currentRead->read->sequence);
            currentRead->reverseRead.reset( new std::string(reverse) );
            int totalKmers = 0;
            currentRead->revKmer.reset( encodeWindow(reverse, &totalKmers) );
            forwardStrand = false;

            SamWriter::Alignment negAlignment;
            std::set<Minhash::Neighbour> negNeighboursCurrentPred = mhIndices[key]->findNeighbours(currentRead->revKmer, *(currentRead->totalKmers));
            bestNeighbourJaccardScore = std::max(bestNeighbourJaccardScore,negNeighboursCurrentPred.begin()->jaccard);
            firstResult = { negNeighboursCurrentPred.begin()->id, forwardStrand};
#if DEBUG_MODE
            LOG(INFO) << "-ve Minhash Neigbours: " << negNeighboursCurrentPred.size() << " in segment: " << std::to_string(pair.first);
#endif
            happy = tryFirstOutofGiven(currentRead, partition, forwardStrand, &referenceGenomeBrigde, negNeighboursCurrentPred, &negAlignment, &score, seenNeighbours);

            if (happy) {
                insertNewAlignment(bestAlignments, {negAlignment, score}, EXTRA_EXTENSION_SHORT);
                seenNeighbours.insert(firstResult.first);
                if (score > bestScore)                 bestScore = score;
                if (bestAlignments.size() >= nBest)    continue;
            }
            queueWrapper.addQueue(pair.first, false, &negNeighboursCurrentPred);

        }
#if DEBUG_MODE
        else {
                    LOG(ERROR) << currentRead->read->key << " predicted as " << std::to_string(pair.first);
                }
#endif
    }

    int counter = 0;
    while(queueWrapper.hasNext() && counter++ < 10 && bestAlignments.size() < nBest) {
        int partition; bool forwardStrand; Minhash::Neighbour neighbour;
        do {
            queueWrapper.pop(&partition, &forwardStrand, &neighbour);
            if (bestNeighbourJaccardScore*0.90 > neighbour.jaccard)
                break;
        } while (queueWrapper.hasNext() && isNeigbourPreviouslySeen(seenNeighbours, neighbour, referenceGenomeBrigde.getWindowLength()));

        bestNeighbourJaccardScore = std::max(bestNeighbourJaccardScore, neighbour.jaccard);

        // if non of the neighbour found by minhash hash atleast 90% score of best terminate the search.
        if (bestNeighbourJaccardScore*0.90 > neighbour.jaccard)
            break;

        SamWriter::Alignment retAlignment; int retScore;
        bool happy = alignMinhashNeighbour(currentRead, partition, forwardStrand, &referenceGenomeBrigde, neighbour,
                                           &retAlignment, &retScore);
        if (happy) {
            insertNewAlignment(bestAlignments, {retAlignment, retScore}, EXTRA_EXTENSION_SHORT);
            seenNeighbours.insert(neighbour.id);
            bestScore = score;
        }
    }
    if ( bestAlignments.size() == 0 ) { // means no mapping found
        SamWriter::Alignment dummAlignment;
        makeUnmapped(currentRead, dummAlignment);
        bestAlignments.push_back({dummAlignment,0});
    }
    return bestAlignments;
}

void align_single(std::string &fastqFile, int &tfBatchSize, TensorflowInference &inferEngine,
                  std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde,
                  SamWriter &samWriter, int nThreads) {
    FastQ fastq(fastqFile);

    while (fastq.hasNext()) {
        std::vector<ReadsWrapper > readsVector(tfBatchSize);
        std::vector< std::pair<std::shared_ptr<Kmer>, int> > pairs(tfBatchSize);
        int loadCount = 0;
        for (int i=0; i<tfBatchSize && fastq.hasNext(); i++) {
            ReadsWrapper readsWrapper;
            readsWrapper.read.reset( new InputRead(fastq.next()) );
            int totalKmers = 0;
            readsWrapper.kmer.reset( encodeWindow((readsWrapper.read)->sequence, &totalKmers) );
            readsWrapper.totalKmers.reset( new int(totalKmers) );
            readsVector[i] = readsWrapper;
            std::pair<std::shared_ptr<Kmer>, int> onePair = { readsVector[i].kmer, totalKmers };
            pairs[i] = onePair;
            loadCount++;
        }

        prediction(inferEngine, readsVector, pairs, loadCount);
        pairs.clear();

        {
            ThreadPool threadPool(nThreads);
            for (int i = 0; i < loadCount; i++) {
                ReadsWrapper eachRead = readsVector[i];
                threadPool.enqueue([eachRead, &mhIndices, &referenceGenomeBrigde, &samWriter] {
                    std::vector<std::pair<SamWriter::Alignment, int>> bestAlignments = alignOne(eachRead, mhIndices, referenceGenomeBrigde, NBEST);
                    for (std::pair<SamWriter::Alignment, int> alig : bestAlignments) {
                        samWriter.writeAlignment(alig.first);
                    }
                });
            }
        }

        LOG(INFO) << loadCount << " finished..";
    }
}


#endif
