//
// Created by Akash Shrestha on 7/16/18.
//

#include <string>
#include "include/tensorflow_inference.h"
#include "include/indexer_jobparser.h"
#include "include/minhash.h"
#include "include/ThreadPool.h"
#include "include/sam_writer.h"
#include "include/input_reader.h"
#include "include/fastaQ.h"
#include "include/paired_priorityqueue_wrapper.h"

struct PairedReadsWrapper {
    std::pair<std::shared_ptr<InputRead>, std::shared_ptr<InputRead> > read;
    std::pair<std::shared_ptr<Kmer>,std::shared_ptr<Kmer> > kmer;
    std::pair<std::shared_ptr<Kmer>,std::shared_ptr<Kmer> > revKmer;
    std::shared_ptr<int> totalKmers;
    std::pair<std::shared_ptr<std::string>, std::shared_ptr<std::string> > reversedRead;
    std::pair<std::shared_ptr<std::set<std::pair<int, bool> > >, std::shared_ptr<std::set<std::pair<int, bool> > > > predictedSegments;

    PairedReadsWrapper() {
        std::nullptr_t null;
        read = {std::shared_ptr<InputRead>(null), std::shared_ptr<InputRead>(null)};
        kmer = { std::shared_ptr<Kmer>(null, [](Kmer* p) { delete [] p; }), std::shared_ptr<Kmer>(null, [](Kmer* p) { delete [] p; }) };
        revKmer = { std::shared_ptr<Kmer>(null, [](Kmer* p) { delete [] p; }), std::shared_ptr<Kmer>(null, [](Kmer* p) { delete [] p; }) };
        totalKmers = std::shared_ptr<int>(null);
        reversedRead = {std::shared_ptr<std::string>(null),std::shared_ptr<std::string>(null)};
        predictedSegments = { std::shared_ptr<std::set<std::pair<int, bool> > >(null), std::shared_ptr<std::set<std::pair<int, bool> > >(null) };
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
    int firstNumMismatches = 0, secondNumMismatches = 0;
    int firstOffset, secondOffset;
    std::string firstPartOfReference = getPartOfReference(refBridge, firstPredictedSegment, firstNeighbour, &firstOffset, windowLength);
    std::string secondPartOfReference = getPartOfReference(refBridge, secondPredictedSegment, secondNeighbour, &secondOffset, windowLength);
    std::string firstQueryString = firstForwardStrand ? currentRead->read.first->sequence : *(currentRead->reversedRead.first);
    std::string secondQueryString = secondForwardStrand ? currentRead->read.second->sequence : *(currentRead->reversedRead.second);

    bool happy;
    SamWriter::Alignment nullAlignment;
    int firstAlignmentScore = SamWriter::alignment(firstPartOfReference, firstQueryString, &nullAlignment, &firstNumMismatches);
    *firstPosition = nullAlignment.pos + firstOffset;

    nullAlignment.pos = 0;
    int secondAlignmentScore = SamWriter::alignment(secondPartOfReference, secondQueryString, &nullAlignment, &secondNumMismatches);
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


inline void pairedPrediction(TensorflowInference &inferEngine, std::vector<PairedReadsWrapper> &readsVector, std::vector< std::pair<std::shared_ptr<Kmer>, int> > &pairs, int loadCount) {
    Tensor tensor = inferEngine.makeTensor(pairs); // pairs.size() = 2*loadCount
    std::vector<std::set<std::pair<int, bool> > > predictions = inferEngine.inference(tensor);
    for (int i=0; i<loadCount; i++) {
        for (int j=0; j<2; j++) {
            if (j%2 == 0) {
                readsVector[i].predictedSegments.first.reset( new std::set<std::pair<int, bool> >() );
                std::set<std::pair<int, bool> > inner = predictions[2 * i + j];
                for (auto lacasito : inner) {
                    readsVector[i].predictedSegments.first->insert(lacasito);
                }
            }
            else {
                readsVector[i].predictedSegments.second.reset( new std::set<std::pair<int, bool> >() );
                std::set<std::pair<int, bool> > inner = predictions[2 * i + j];
                for (auto lacasito : inner) {
                    readsVector[i].predictedSegments.second->insert(lacasito);
                }
            }
        }
    }
    predictions.clear();
}

std::string alignOne_paired(PairedReadsWrapper eachRead, std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde) {
    PairedPriorityQueueWrapper queueWrapper;
    PairedReadsWrapper *currentRead = &eachRead;
    uint32_t firstReadPosition = 0, secondReadPosition = 0;
    bool firstAlignmentForward = false, secondAlignmentForward = false;
    int bestScore = -1 * (currentRead->read.first->sequence.length()+currentRead->read.second->sequence.length());
    int score = bestScore;
    for (std::pair<int, bool> firstPrediction : *(currentRead->predictedSegments.first)) {
        if (firstPrediction.first >= mhIndices.size())
            continue;
        std::string key1 = "index-" + std::to_string(firstPrediction.first) + ".mh";
        std::set<Minhash::Neighbour> firstPosNeighbours = mhIndices[key1]->findNeighbours(currentRead->kmer.first, *(currentRead->totalKmers));

        string firstReverse = reverseComplement(currentRead->read.first->sequence);
        currentRead->reversedRead.first.reset( new std::string(firstReverse) );
        int totalKmers = 0;
        currentRead->revKmer.first.reset( encodeWindow(firstReverse, &totalKmers) );
        std::set<Minhash::Neighbour> firstNegNeighbours = mhIndices[key1]->findNeighbours(currentRead->revKmer.first, *(currentRead->totalKmers));

        for (std::pair<int, bool> secondPrediction : *(currentRead->predictedSegments.second)) {
            if (secondPrediction.first >= mhIndices.size() || abs(firstPrediction.first-secondPrediction.first)>1 )
                continue;
            std::string key2 = "index-" + std::to_string(secondPrediction.first) + ".mh";
            std::set<Minhash::Neighbour> secondPosNeighbours = mhIndices[key2]->findNeighbours(currentRead->kmer.second, *(currentRead->totalKmers));

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

            std::string secondReverse = reverseComplement(currentRead->read.second->sequence);
            currentRead->reversedRead.second.reset( new std::string(secondReverse) );
            int totalKmers = 0;
            currentRead->revKmer.second.reset( encodeWindow(secondReverse, &totalKmers) );
            std::set<Minhash::Neighbour> secondNegNeighbours = mhIndices[key2]->findNeighbours(currentRead->revKmer.second, *(currentRead->totalKmers));
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

    std::stringstream stringstream;
    stringstream << currentRead->read.first->key << std::endl;
    stringstream << "firstPosition: " << firstReadPosition << " forward:" << firstAlignmentForward
              << " secondPosition: " << secondReadPosition << " forward: " << secondAlignmentForward << std::endl;
    return stringstream.str();

}

std::mutex mutex;
void dealWithAlignment(std::string &alignment) {
    std::unique_lock<std::mutex> lock (mutex);
    std::cout << alignment;
}

void align_paired(std::string &fastqFiles, int &tfBatchSize, TensorflowInference &inferEngine,
                  std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde,
                  SamWriter &samWriter, int nThreads) {
    std::vector<std::string> files = split(fastqFiles, ",");
    if (files.size() != 2) {
        throw new InputReaderException("paired alignment files should be \"<1stFile>,<2ndFile>\"");
    }

    int K = referenceGenomeBrigde.getK();
    std::string fastQFiles[2] = {files[0], files[1]};
    PairedFastQ fastQ(fastQFiles);
    while (fastQ.hasNext()) {
        std::vector<PairedReadsWrapper> readsVector(tfBatchSize);
        std::vector< std::pair<std::shared_ptr<Kmer>, int> > pairs(tfBatchSize*2);
        int loadCount = 0;
        for (int i=0; i<tfBatchSize && fastQ.hasNext(); i++) {
            PairedReadsWrapper readsWrapper;
            InputRead pairedReads[2];
            fastQ.next(pairedReads);
            readsWrapper.read.first.reset( new InputRead(pairedReads[0]) );
            readsWrapper.read.second.reset( new InputRead(pairedReads[1]) );
            int totalKmer = 0;
            readsWrapper.kmer.first.reset( encodeWindow(readsWrapper.read.first->sequence, &totalKmer, K) );
            readsWrapper.kmer.second.reset( encodeWindow(readsWrapper.read.second->sequence, &totalKmer, K) );
            readsWrapper.totalKmers.reset(new int(totalKmer));
            readsVector[i] = readsWrapper;
            pairs[2*i] = std::pair<std::shared_ptr<Kmer>, int> (readsWrapper.kmer.first, totalKmer);
            pairs[2*i+1] = std::pair<std::shared_ptr<Kmer>, int> (readsWrapper.kmer.second, totalKmer);
            loadCount++;
        }

        pairedPrediction(inferEngine, readsVector, pairs, loadCount);
        pairs.clear();

        {
            ThreadPool threadPool(nThreads);
            for (int i = 0; i < loadCount; i++) {
                PairedReadsWrapper eachRead = readsVector[i];
                threadPool.enqueue([eachRead, &mhIndices, &referenceGenomeBrigde, &samWriter] {
                    std::string alignment = alignOne_paired(eachRead, mhIndices, referenceGenomeBrigde);
                    dealWithAlignment(alignment);
                });
            }
        }
        std::cout << loadCount << " done" << std::endl;
    }
}