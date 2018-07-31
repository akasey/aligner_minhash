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
    std::auto_ptr<InputRead> read[2];
    Kmer * kmer[2] = {NULL, NULL};
    Kmer * revKmer[2] = {NULL, NULL};
    std::auto_ptr<int> totalKmer;
    std::auto_ptr<std::string> reversedRead[2];
    std::auto_ptr<std::set<std::pair<int, bool> > >predictedSegments[2];

    void clear() {
        for (int j=0; j<2; j++) {
            if (kmer[j]) delete[] kmer[j];
            if (revKmer[j]) delete[] revKmer[j];
        }
    }
};

inline std::string getPartOfReference(IndexerJobParser *refBridge, int &predictedSegment,
                                      Minhash::Neighbour &neighbour, int *startingCoord, int &windowLength) {
    std::pair<int, std::string> *refSegment = refBridge->getSegmentForID(predictedSegment);
    NULL_CHECK(refSegment, "Reference for segment " + std::to_string(predictedSegment) + " is NULL");
    int extraLength = (int) ((float)windowLength * 0.10);
    int start = fmax((int)(neighbour.id) - (int)(refSegment->first) - extraLength, 0);
    int length = fmin(refSegment->second.length(), start+2*extraLength+windowLength) - start;
    std::string partOfReference = refSegment->second.substr(start, length);
    *startingCoord = start + refSegment->first;
    return partOfReference;
}

inline bool alignPairedMinhashNeighbour(PairedReadsWrapper *currentRead, int &firstPredictedSegment, bool &firstForwardStrand,
            int &secondPredictedSegment, bool &secondForwardStrand, IndexerJobParser *refBridge,
            Minhash::Neighbour &firstNeighbour, Minhash::Neighbour &secondNeighbour,
            int windowLength, int *score, uint32_t *firstPosition, uint32_t *secondPosition) {
    int firstOffset, secondOffset;
    std::string firstPartOfReference = getPartOfReference(refBridge, firstPredictedSegment, firstNeighbour, &firstOffset, windowLength);
    std::string secondPartOfReference = getPartOfReference(refBridge, secondPredictedSegment, secondNeighbour, &secondOffset, windowLength);
    std::string firstQueryString = firstForwardStrand ? currentRead->read[0]->sequence : *(currentRead->reversedRead[0]);
    std::string secondQueryString = secondForwardStrand ? currentRead->read[1]->sequence : *(currentRead->reversedRead[1]);

    bool happy;
    SamWriter::Alignment nullAlignment;
    int firstAlignmentScore = SamWriter::alignment(firstPartOfReference, firstQueryString, &nullAlignment);
    *firstPosition = nullAlignment.pos + firstOffset;

    nullAlignment.pos = 0;
    int secondAlignmentScore = SamWriter::alignment(secondPartOfReference, secondQueryString, &nullAlignment);
    *secondPosition = nullAlignment.pos + secondOffset;
    *score = (firstAlignmentScore+secondAlignmentScore);

    if ( (firstQueryString.length()+secondQueryString.length())*0.8 <= *score ) {
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


inline void pairedPrediction(TensorflowInference &inferEngine, std::vector<PairedReadsWrapper> &readsVector, std::vector< std::pair<Kmer *, int> > &pairs, int loadCount) {
    Tensor tensor = inferEngine.makeTensor(pairs); // pairs.size() = 2*loadCount
    std::vector<std::set<std::pair<int, bool> > > predictions = inferEngine.inference(tensor);
    for (int i=0; i<loadCount; i++) {
        for (int j=0; j<2; j++) {
            readsVector[i].predictedSegments[j] = std::auto_ptr<std::set<std::pair<int, bool> >>(new std::set<std::pair<int, bool> > );
            std::set<std::pair<int, bool> > inner = predictions[2*i+j];
            for (auto lacasito : inner) {
                readsVector[i].predictedSegments[j]->insert(lacasito);
            }
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

    int K = referenceGenomeBrigde.getK();
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
            readsWrapper.read[0] = std::auto_ptr<InputRead>(new InputRead(pairedReads[0]));
            readsWrapper.read[1] = std::auto_ptr<InputRead>(new InputRead(pairedReads[1]));
            int totalKmer = 0;
            readsWrapper.kmer[0] = encodeWindow(readsWrapper.read[0]->sequence, &totalKmer, K);
            readsWrapper.kmer[1] = encodeWindow(readsWrapper.read[1]->sequence, &totalKmer, K);
            readsWrapper.totalKmer = std::auto_ptr<int>(new int(totalKmer));
            readsVector[i] = readsWrapper;
            std::pair<Kmer *, int> onePair(readsWrapper.kmer[0], totalKmer);
            std::pair<Kmer *, int> twoPair(readsWrapper.kmer[1], totalKmer);
            pairs[2*i] = onePair;
            pairs[2*i+1] = twoPair;
            loadCount++;
        }

        pairedPrediction(inferEngine, readsVector, pairs, loadCount);

#if 0
        for (int i=0; i<loadCount; i++) {
            PairedReadsWrapper *currentRead = &(readsVector[i]);
            for (int j=0; j<2; j++) {
                std::cout << j << "th ";
                for (auto lacasito : *(currentRead->predictedSegments[j])) {
                    std::cout << lacasito.first << " ";
                }
            }
            std::cout << " \n";
        }
        continue;
#endif

        PairedPriorityQueueWrapper queueWrapper;
        for (int i=0; i<loadCount; i++) {
            PairedReadsWrapper *currentRead = &(readsVector[i]);
            uint32_t firstReadPosition = 0, secondReadPosition = 0;
            bool firstAlignmentForward = false, secondAlignmentForward = false;
            int bestScore = -1 * (currentRead->read[0]->sequence.length()+currentRead->read[1]->sequence.length());
            int score = bestScore;
            for (std::pair<int, bool> firstPrediction : *(currentRead->predictedSegments[0])) {
                if (firstPrediction.first >= mhIndices.size())
                    continue;
                std::string key1 = "index-" + std::to_string(firstPrediction.first) + ".mh";
                std::set<Minhash::Neighbour> firstPosNeighbours = mhIndices[key1]->findNeighbours(currentRead->kmer[0], *(currentRead->totalKmer));

                string firstReverse = reverseComplement(currentRead->read[0]->sequence);
                currentRead->reversedRead[0] = std::auto_ptr<std::string>(new std::string(firstReverse));
                int totalKmers = 0;
                currentRead->revKmer[0] = encodeWindow(firstReverse, &totalKmers);
                std::set<Minhash::Neighbour> firstNegNeighbours = mhIndices[key1]->findNeighbours(currentRead->revKmer[0], *(currentRead->totalKmer));

                for (std::pair<int, bool> secondPrediction : *(currentRead->predictedSegments[1])) {
                    if (secondPrediction.first >= mhIndices.size() || abs(firstPrediction.first-secondPrediction.first)>1 )
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
                            firstAlignmentForward = firstStrandForward; secondAlignmentForward = secondStrandForward;
                            bestScore = score;
                            continue;
                        }
                    }

                    std::string secondReverse = reverseComplement(currentRead->read[1]->sequence);
                    currentRead->reversedRead[1] = std::auto_ptr<std::string>(new std::string(secondReverse));
                    int totalKmers = 0;
                    currentRead->revKmer[1] = encodeWindow(secondReverse, &totalKmers);
                    std::set<Minhash::Neighbour> secondNegNeighbours = mhIndices[key2]->findNeighbours(currentRead->revKmer[1], *(currentRead->totalKmer));
                    firstStrandForward = true; secondStrandForward = false;
                    happy = tryPairedFirstOutOfGiven(currentRead, firstPredictedSegment, firstStrandForward,
                                                          secondPredictedSegment, secondStrandForward, &referenceGenomeBrigde,
                                                          firstPosNeighbours, secondNegNeighbours,
                                                          referenceGenomeBrigde.getWindowLength(), &score, &firstPosition, &secondPosition);
                    if (!happy) {
                        queueWrapper.addQueue(firstPrediction.first, firstStrandForward, &firstPosNeighbours, secondPrediction.first, secondStrandForward, &secondNegNeighbours);
                    }
                    else {
                        if (score > bestScore) {
                            firstReadPosition = firstPosition;
                            secondReadPosition = secondPosition;
                            firstAlignmentForward = firstStrandForward; secondAlignmentForward = secondStrandForward;
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

        for (int i=0; i<loadCount; i++) {
            readsVector[i].clear();
        }
        std::cout << loadCount << " done" << std::endl;
    }
}