//
// Created by Akash Shrestha on 7/16/18.
//

#include <string>
#include "include/tensorflow_inference.h"
#include "include/indexer_jobparser.h"
#include "include/minhash.h"
#include "include/sam_writer.h"
#include "include/input_reader.h"
#include "include/fastaQ.h"
#include "include/paired_priorityqueue_wrapper.h"

struct PairedReadsWrapper {
    InputRead *read[2];
    Kmer *kmer[2];
    Kmer *revKmer[2];
    int *totalKmer;
    std::string *reversedRead[2];
    std::set<std::pair<int, bool>> *predictedSegments[2];

    void clear() {
        if (totalKmer) delete totalKmer;
        for (int i=0; i<2; i++) {
            if (read[i]) delete[] read[i];
            if (kmer[i]) delete[] kmer[i];
            if (revKmer[i]) delete[] revKmer[i];
            if (reversedRead[i]) delete reversedRead[i];
            if (predictedSegments[i]) delete predictedSegments[i];
        }
    }
};

inline std::string getPartOfReference(IndexerJobParser *refBridge, int &predictedSegment,
                                      Minhash::Neighbour &neighbour, int &windowLength) {
    std::pair<int, std::string> *refSegment = refBridge->getSegmentForID(predictedSegment);
    NULL_CHECK(refSegment, "Reference for segment " + std::to_string(predictedSegment) + " is NULL");
    int extraLength = (int) ((float)windowLength * 0.10);
    int start = fmax((int)(neighbour.id) - (int)(refSegment->first) - extraLength, 0);
    int length = fmin(refSegment->second.length(), start+2*extraLength+windowLength) - start;
    std::string partOfReference = refSegment->second.substr(start, length);
    return partOfReference;
}

inline bool alignPairedMinhashNeighbour(PairedReadsWrapper *currentRead, int &firstPredictedSegment, bool &firstForwardStrand,
            int &secondPredictedSegment, bool &secondForwardStrand, IndexerJobParser *refBridge,
            Minhash::Neighbour &firstNeighbour, Minhash::Neighbour &secondNeighbour,
            int windowLength, int *score, uint32_t *firstPosition, uint32_t *secondPosition) {
    std::string firstPartOfReference = getPartOfReference(refBridge, firstPredictedSegment, firstNeighbour, windowLength);
    std::string secondPartOfReference = getPartOfReference(refBridge, secondPredictedSegment, secondNeighbour, windowLength);
    std::string firstQueryString = firstForwardStrand ? currentRead->read[0]->sequence : *(currentRead->reversedRead[0]);
    std::string secondQueryString = secondForwardStrand ? currentRead->read[1]->sequence : *(currentRead->reversedRead[1]);

    bool happy;
    SamWriter::Alignment nullAlignment;
    int firstAlignmentScore = SamWriter::alignment(firstPartOfReference, firstQueryString, &nullAlignment);
    *firstPosition = nullAlignment.pos;

    nullAlignment.pos = 0;
    int secondAlignmentScore = SamWriter::alignment(secondPartOfReference, secondQueryString, &nullAlignment);
    *firstPosition = nullAlignment.pos;

    if ( (firstQueryString.length()+secondQueryString.length())*0.6 <= (firstAlignmentScore+secondAlignmentScore) ) {
        happy = true;
    }
    else {
        happy = false;
    }
    return happy;
}

inline bool tryPairedFirstOutOfGiven(PairedReadsWrapper *currentRead, int &firstPredictedSegment, bool &firstForwardStrand,
                               int &secondPredictedSegment, bool &secondForwardStrand, IndexerJobParser *refBridge,
                               std::set<Minhash::Neighbour> &firstNeighbour, std::set<Minhash::Neighbour> &secondNeighbour,
                               int windowLength, int *score, uint32_t *firstPosition, uint32_t *secondPosition) {
    if (firstNeighbour.size() > 0 && secondNeighbour.size() > 0) {
        Minhash::Neighbour first = *(firstNeighbour.begin());
        firstNeighbour.erase(firstNeighbour.begin());
        Minhash::Neighbour second = *(secondNeighbour.begin());
        secondNeighbour.erase(secondNeighbour.begin());
        return alignPairedMinhashNeighbour(currentRead, firstPredictedSegment, firstForwardStrand,
                        secondPredictedSegment, secondForwardStrand, refBridge, first, second,
                                           windowLength, score, firstPosition, secondPosition);
    }
    return false;
}


inline void pairedPrediction(TensorflowInference &inferEngine, std::vector<PairedReadsWrapper> &readsVector, std::vector< std::pair<Kmer*,int> > &pairs, int loadCount) {
    Tensor tensor = inferEngine.makeTensor(pairs); // pairs.size() = 2*loadCount
    std::vector<std::set<std::pair<int, bool> > > predictions = inferEngine.inference(tensor);
    for (int i=0; i<loadCount; i++) {
        readsVector[i].predictedSegments[0] = new std::set<std::pair<int, bool> >();
        readsVector[i].predictedSegments[1] = new std::set<std::pair<int, bool> >();
        std::set<std::pair<int, bool> > inner = predictions[i];
        int j = 0;
        for (auto lacasito : inner) {
            int key = j%2;
            readsVector[i].predictedSegments[key]->insert(lacasito);
            j++;
        }
    }
    predictions.clear();
}

void align_paired(std::string &fastqFiles, int &tfBatchSize, TensorflowInference &inferEngine,
                  std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde,
                  SamWriter &samWriter) {
    std::vector<std::string> files = split(fastqFiles, ",");
    if (files.size() != 2) {
        throw new InputReaderException("paired alignment files should be \"<1stFile>,<2ndFile>\"");
    }

    std::string fastQFiles[2] = {files[0], files[1]};
    PairedFastQ fastQ(fastQFiles);
    while (fastQ.hasNext()) {
        std::vector<PairedReadsWrapper> readsVector(tfBatchSize);
        std::vector< std::pair<Kmer*, int> > pairs(tfBatchSize*2);
        int loadCount = 0;
        for (int i=0; i<tfBatchSize && fastQ.hasNext(); i++) {
            PairedReadsWrapper readsWrapper;
            InputRead pairedReads[2];
            fastQ.next(pairedReads);
            readsWrapper.read[0] = new InputRead(pairedReads[0]);
            readsWrapper.read[1] = new InputRead(pairedReads[1]);
            int totalKmer = 0;
            readsWrapper.kmer[0] = encodeWindow(readsWrapper.read[0]->sequence, &totalKmer);
            readsWrapper.kmer[1] = encodeWindow(readsWrapper.read[1]->sequence, &totalKmer);
            readsWrapper.totalKmer = new int(totalKmer);
            readsVector[i] = readsWrapper;
            std::pair<Kmer *, int> onePair(readsVector[i].kmer[0], totalKmer);
            std::pair<Kmer *, int> twoPair(readsVector[i].kmer[1], totalKmer);
            pairs[2*i] = onePair;
            pairs[2*i+1] = twoPair;
            loadCount++;
        }

        pairedPrediction(inferEngine, readsVector, pairs, loadCount);

        PairedPriorityQueueWrapper queueWrapper;
        for (int i=0; i<loadCount; i++) {
            PairedReadsWrapper *currentRead = &(readsVector[i]);
            uint32_t firstReadPosition, secondReadPosition;
            bool firstAlignmentForward, secondAlignmentForward;
            int bestScore = -1 * (currentRead->read[0]->sequence.length()+currentRead->read[1]->sequence.length());
            int score;
            for (std::pair<int, bool> firstPrediction : *(currentRead->predictedSegments[0])) {
                if (firstPrediction.first >= mhIndices.size())
                    continue;
                std::string key1 = "index-" + std::to_string(firstPrediction.first) + ".mh";
                std::set<Minhash::Neighbour> firstPosNeighbours = mhIndices[key1]->findNeighbours(currentRead->kmer[0], *(currentRead->totalKmer));

                string firstReverse = reverseComplement(currentRead->read[0]->sequence);
                currentRead->reversedRead[0] = new std::string(firstReverse);
                int totalKmers = 0;
                currentRead->revKmer[0] = encodeWindow(firstReverse, &totalKmers);
                std::set<Minhash::Neighbour> firstNegNeighbours = mhIndices[key1]->findNeighbours(currentRead->revKmer[0], *(currentRead->totalKmer));

                for (std::pair<int, bool> secondPrediction : *(currentRead->predictedSegments[1])) {
                    if (secondPrediction.first >= mhIndices.size() )//|| abs(firstPrediction.first-secondPrediction.first)>1 )
                        continue;
                    std::string key2 = "index-" + std::to_string(secondPrediction.first) + ".mh";
                    std::set<Minhash::Neighbour> secondPosNeighbours = mhIndices[key2]->findNeighbours(currentRead->kmer[1], *(currentRead->totalKmer));

                    int firstPredictedSegment = firstPrediction.first, secondPredictedSegment = secondPrediction.first;
                    uint32_t firstPosition, secondPosition;
                    bool firstStrandForward = false, secondStrandForward = true;

                    bool happy = tryPairedFirstOutOfGiven(currentRead, firstPredictedSegment, firstStrandForward,
                                                          secondPredictedSegment, secondStrandForward, &referenceGenomeBrigde,
                                                          firstNegNeighbours, secondPosNeighbours,
                                    referenceGenomeBrigde.getWindowLength(), &score, &firstPosition, &secondPosition);
                    if (!happy) {
                        queueWrapper.addQueue(firstPrediction.first, firstStrandForward, &firstNegNeighbours, secondPrediction.first, secondStrandForward, &secondPosNeighbours);
                    }
                    else {
                        if (score > bestScore) {
                            firstReadPosition = firstPosition;
                            secondReadPosition = secondPosition;
                            firstAlignmentForward = false; secondAlignmentForward = true;
                            bestScore = score;
                            continue;
                        }
                    }

                    std::string secondReverse = reverseComplement(currentRead->read[1]->sequence);
                    currentRead->reversedRead[1] = new std::string(secondReverse);
                    int totalKmers = 0;
                    currentRead->revKmer[1] = encodeWindow(secondReverse, &totalKmers);
                    std::set<Minhash::Neighbour> secondNegNeighbours = mhIndices[key2]->findNeighbours(currentRead->revKmer[1], *(currentRead->totalKmer));
                    happy = tryPairedFirstOutOfGiven(currentRead, firstPredictedSegment, firstStrandForward,
                                                          secondPredictedSegment, secondStrandForward, &referenceGenomeBrigde,
                                                          firstNegNeighbours, secondPosNeighbours,
                                                          referenceGenomeBrigde.getWindowLength(), &score, &firstPosition, &secondPosition);
                    firstStrandForward = true; secondStrandForward = false;
                    if (!happy) {
                        queueWrapper.addQueue(firstPrediction.first, firstStrandForward, &firstPosNeighbours, secondPrediction.first, secondStrandForward, &secondNegNeighbours);
                    }
                    else {
                        if (score > bestScore) {
                            firstReadPosition = firstPosition;
                            secondReadPosition = secondPosition;
                            firstAlignmentForward = true; secondAlignmentForward = false;
                            bestScore = score;
                            continue;
                        }
                    }
                }
            }

            while (queueWrapper.hasNext()) {
                uint32_t firstPosition, secondPosition;
                int firstPartition; bool firstForwardStrand; Minhash::Neighbour firstNeighbour;
                int secondPartition; bool secondForwardStrand; Minhash::Neighbour secondNeighbour;
                queueWrapper.pop(&firstPartition, &firstForwardStrand, &firstNeighbour,
                                 &secondPartition, &secondForwardStrand, &secondNeighbour);
                bool happy = alignPairedMinhashNeighbour(currentRead, firstPartition, firstForwardStrand,
                                                 secondPartition, secondForwardStrand, &referenceGenomeBrigde,
                                                 firstNeighbour, secondNeighbour,
                                                 referenceGenomeBrigde.getWindowLength(), &score, &firstPosition, &secondPosition);
                if (happy) {
                    firstReadPosition = firstPosition; secondReadPosition = secondPosition;
                    firstAlignmentForward = firstForwardStrand; secondAlignmentForward = secondForwardStrand;
                    break;
                }
            }

            std::cout << currentRead->read[0]->key << std::endl;
            std::cout << "firstPosition: " << firstReadPosition << " forward:" << firstAlignmentForward
                      << " secondPosition: " << secondReadPosition << " forward: " << secondAlignmentForward << std::endl;

        }
        std::cout << loadCount << " done" << std::endl;
    }
}