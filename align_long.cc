//
// Created by Akash Shrestha on 9/12/18.
//

/*
 * Single read aligner for third generation reads which are longer than NGS.
 * Break then into shorter numFragments and predict
 */

#ifndef ALIGNER_ALIGN_LONG_H
#define ALIGNER_ALIGN_LONG_H

#include "include/fixed_array.h"
#include "include/sam_writer.h"
#include "include/common.h"
#include "include/minhash.h"
#include "include/tf_metaparser.h"
#include "include/tensorflow_inference.h"
#include "include/ThreadPool.h"
#include "include/fastaQ.h"
#include "include/sam_writer.h"
#include "include/indexer_jobparser.h"
#include "include/long_priorityqueue_wrapper.h"

#include <list>


struct Prediction {
    std::vector<std::pair<int, int> > predictions; // (fragment, predictedSegment)
    bool forward = true;

    void insert(int fragment, int predictedSegment) {
        predictions.push_back({fragment, predictedSegment});
    }
    bool operator<(const Prediction &rhs) const {
        return predictions.size() > rhs.predictions.size();
    }
    void print(std::ostream &out) {
        out << "On-forward:" << (forward ? "true":"false") << " ";
        for (const auto &each : predictions) {
            out << each.first << ":" << each.second << " ";
        }
        out << std::endl;
    }
};

struct LongReadsWrapper {
    std::shared_ptr<InputRead> read;
    std::shared_ptr<std::string> reverseRead;
    std::shared_ptr<int> numFragments;
    std::shared_ptr<std::vector<std::string> > positiveFragments;
    std::shared_ptr<std::vector<std::shared_ptr<Kmer> > > kmers;
    std::shared_ptr<std::vector<std::shared_ptr<Kmer> > > revKmers;
    std::shared_ptr<int> totalKmers;
    std::shared_ptr<std::vector<Prediction > > predictions;

    LongReadsWrapper() {
    }

};

inline std::string extension(int direction, IndexerJobParser *refBridge, int segment, int startPoint, int requiredLength, std::list<std::tuple<int,int,int> > *pathOfSegments) {
    try {
        std::pair<int, std::string> *referenceSegment = refBridge->getSegmentForID(segment);
        if (referenceSegment == NULL || requiredLength == 0) return "";
        std::string substr = "";
        if (direction > 0) { // right extension
            int startCoord = startPoint - referenceSegment->first;
            if (startCoord < 0) return "";
            int length = fmin(referenceSegment->second.length() - startCoord, requiredLength);
            substr = referenceSegment->second.substr(startCoord, length);
            pathOfSegments->push_back(std::make_tuple(segment, startPoint, startPoint + length));
            if (requiredLength - length > 0) {
                substr = substr + extension(direction, refBridge, segment + 1, startPoint + length, requiredLength - length, pathOfSegments);
            }
        } else { // left extension
            int startCoord = startPoint - referenceSegment->first;
            if (startCoord < 0) return "";
            int start = fmax(startCoord - requiredLength, 0);
            int length = startCoord - start;
            substr = referenceSegment->second.substr(start, length);
            pathOfSegments->push_front(std::make_tuple(segment, referenceSegment->first+start, referenceSegment->first+start+length));
            if (requiredLength - length > 0) {
                substr = extension(direction, refBridge, segment - 1, startPoint - length, requiredLength - length, pathOfSegments) + substr;
            }
        }
        return substr;
    }
    catch (std::out_of_range &ex) {
        return "";
    }
}

inline std::string getReadFromReference(int fragment, int segment, int length, int startPoint, bool forwardStrand, IndexerJobParser *refBridge, std::list<std::tuple<int,int,int> > *pathOfSegments) {
    int startPositionForFragment = (refBridge->getWindowLength()/2)*fragment;
    int requiredLengthAtLeft = startPositionForFragment + 500;
    int requiredLengthAtRight = length-startPositionForFragment + 500;
    if (!forwardStrand) {
        int temp = requiredLengthAtLeft;
        requiredLengthAtLeft = requiredLengthAtRight;
        requiredLengthAtRight = temp;
    }
    std::string referenceSegment = extension(-1, refBridge, segment, startPoint, requiredLengthAtLeft, pathOfSegments)
            + extension(1, refBridge, segment, startPoint, requiredLengthAtRight, pathOfSegments);
    return referenceSegment;
}

inline bool alignMinhashNeighbour_long(LongReadsWrapper *currentRead, int &fragment,
                                  int &predictedSegment, bool &forwardStrand, IndexerJobParser *refBridge, Minhash::Neighbour &neighbour,
                                  SamWriter::Alignment *retAlignment, int *score) {
    std::list<std::tuple<int,int,int> > pathOfSegment; // segment, start, end
    std::string referenceSegment = getReadFromReference(fragment, predictedSegment, currentRead->read->sequence.length(), neighbour.id, forwardStrand, refBridge, &pathOfSegment);
    int start = std::get<1>(pathOfSegment.front());
    int numMismatches = 0;

    std::string queryString = forwardStrand ? currentRead->read->sequence : *(currentRead->reverseRead);
    *score = SamWriter::alignment(referenceSegment, queryString, retAlignment, &numMismatches);

#if DEBUG_MODE
    std::cout << "Reference: " << std::endl;
    std::cout << referenceSegment << std::endl;
    std::cout << "QueryString: " << std::endl;
    std::cout << queryString << std::endl;
#endif

    bool happy;
    retAlignment->qname = currentRead->read->key;
    std::string segmentName = refBridge->segIdToKey(predictedSegment);
//    if (segmentName.empty()) throw
    retAlignment->rname = split(segmentName, " ")[0];
    retAlignment->pos = retAlignment->pos + start-1;

    // Can't rely on score of alignment. Need to see alignment cigar
//    std::retAlignment->cigar
#if DEBUG_MODE
    LOG(INFO) << "Alignment score: " << *score << " Happy threshold: " << queryString.length()*0.65;
#endif
//    if ( queryString.length()*0.65 <= *score) { // atleast 70% matches // consider mapped
    if ( numMismatches <= 0.10 * queryString.length() ) { // atleast 70% matches // consider mapped
        happy = true;
        retAlignment->flag = retAlignment->flag | (forwardStrand ? 0 : REVERSE_MAPPED);
    }
    else {
        happy = false;
        retAlignment->flag = retAlignment->flag | SEGMENT_UNMAPPED;
    }
#if DEBUG_MODE
    retAlignment->print(std::cout);
    std::cout << *score << "/" << currentRead->read->sequence.length() << std::endl;
#endif
    return happy;
}

inline bool tryFirstNeighbour(LongReadsWrapper *currentRead, int &fragment,
                              int &predictedSegment, bool &forwardStrand, IndexerJobParser *refBridge, std::set<Minhash::Neighbour> &neighbours,
                              SamWriter::Alignment *retAlignment, int *score) {
    if (neighbours.size() > 0) {
        Minhash::Neighbour first = *(neighbours.begin());
        neighbours.erase(neighbours.begin());
        return alignMinhashNeighbour_long(currentRead, fragment, predictedSegment, forwardStrand, refBridge, first, retAlignment, score);
    }
    return false;
}

inline void prediction_long(TensorflowInference &inferEngine, std::vector<LongReadsWrapper > &readsVector, int totalFragments, int loadCount, size_t mhIndicesSize) {
    std::vector< std::pair<std::shared_ptr<Kmer>, int> > pairs(totalFragments);
    int filler = 0;
    for (int i=0; i<loadCount; i++) {
        LongReadsWrapper *currentRead = &readsVector[i];
        for (int j=0; j<*(currentRead->numFragments.get()); j++) {
            pairs[filler++] = {currentRead->kmers->at(j), *(currentRead->totalKmers)};
        }
    }

    Tensor tensor = inferEngine.makeTensor(pairs);
    std::vector<std::set<std::pair<int, bool> > > predictions = inferEngine.inference(tensor);
    pairs.clear();

    filler = 0;
    for (int i=0; i<loadCount; i++) {
        LongReadsWrapper *currentRead = &readsVector[i];
        int numFragments = *(currentRead->numFragments);
        std::vector<std::set<std::pair<int,bool> > > thisLoadPrediction(numFragments);
        std::vector<Prediction> *predictSequence = new std::vector<Prediction>();
        for (int j=0; j<numFragments; j++) {
            std::set<std::pair<int, bool> > copy;
            for (std::pair<int,bool> each  : predictions[filler]) {
                if (each.first < mhIndicesSize) {
                    copy.insert(each);
                }
            }
            thisLoadPrediction.at(j) = copy;
            filler++;
        }
        bool noMore = false;
        while (!noMore) {
            noMore = true;
            Prediction pred;
            int lastSequence = -1;
            for (int j=0; j<numFragments; j++) {
                if (thisLoadPrediction[j].size() > 0) {
                    std::set<std::pair<int,bool> >::iterator top = thisLoadPrediction[j].begin();
                    if ( lastSequence == -1 || abs(lastSequence - top->first) <=1 ) {
                        pred.insert(j, top->first);
                        if (lastSequence - top->first > 0)
                            pred.forward = false;
                        lastSequence = top->first;
                        thisLoadPrediction[j].erase(top);
                    }
                }
            }
            if (pred.predictions.size() != 0) {
                predictSequence->push_back(pred);
                noMore = false;
            }

        }
        std::sort(predictSequence->begin(), predictSequence->end());
        currentRead->predictions.reset(predictSequence);
    }
}

inline void makeUnmapped_long(LongReadsWrapper *currentRead, SamWriter::Alignment &retAlignment) {
    retAlignment.qname = currentRead->read->key;
    retAlignment.flag = retAlignment.flag | SEGMENT_UNMAPPED;
}

SamWriter::Alignment alignOne_long(LongReadsWrapper eachRead, std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde) {
    LongPriorityQueueWrapper queueWrapper;
    LongReadsWrapper *currentRead = &eachRead;
    SamWriter::Alignment bestAlignment;
    int bestScore = -1 * (currentRead->read->sequence.length());
    int score = 0;
    bool foundOne = false;
    for (int j=0; j<currentRead->predictions->size() && !foundOne; j++) {
#if DEBUG_MODE
        std::cout << "Predictions " << std::endl;
        currentRead->predictions->at(j).print(std::cout);
#endif
        Prediction *currentPrediction = &(currentRead->predictions->at(j));
        for (const std::pair<int, int> &fragSeg : currentPrediction->predictions) {
            std::string key = "index-" + std::to_string(fragSeg.second) + ".mh";

            // Schedule positive to go first or negative depending on prediction sequence is decreasing or increasing
            std::shared_ptr<std::vector<std::shared_ptr<Kmer> > > order1, order2;
            bool forward1, forward2;
            if (currentPrediction->forward) {
                order1 = currentRead->kmers; order2 = currentRead->revKmers;
                forward1 = true; forward2 = false;
            }
            else {
                order1 = currentRead->revKmers; order2 = currentRead->kmers;
                forward1 = false; forward2 = true;
            }

            std::set<Minhash::Neighbour> order1NeighboursCurrentPred = mhIndices[key]->findNeighbours(order1->at(fragSeg.first), *(currentRead->totalKmers.get()));
            SamWriter::Alignment alignment;
            int fragment = fragSeg.first, segment = fragSeg.second;
            bool forwardStrand = forward1;
            bool happy = tryFirstNeighbour(currentRead, fragment, segment, forwardStrand, &referenceGenomeBrigde,
                                           order1NeighboursCurrentPred, &alignment, &score);
            if (!happy) {
                queueWrapper.addQueue(fragment, segment, forwardStrand, &order1NeighboursCurrentPred);
            }
            else if (score > bestScore){
                bestAlignment = alignment;
                bestScore = score;
                foundOne = true;
                break;
            }

            // Try negative strand
            std::set<Minhash::Neighbour> order2NeighboursCurrentPred = mhIndices[key]->findNeighbours(order2->at(fragSeg.first), *(currentRead->totalKmers.get()));
            SamWriter::Alignment negAlignment;
            forwardStrand = forward2;
            happy = tryFirstNeighbour(currentRead, fragment, segment, forwardStrand, &referenceGenomeBrigde,
                                      order2NeighboursCurrentPred, &negAlignment, &score);
            if (!happy) {
                queueWrapper.addQueue(fragment, segment, forwardStrand, &order2NeighboursCurrentPred);
            }
            else if (score > bestScore){
                bestAlignment = negAlignment;
                bestScore = score;
                foundOne = true;
                break;
            }
        }
    }

    int counter = 0;
    while (queueWrapper.hasNext() && !foundOne && counter++ < 10) {
        int fragment, partition; bool forwardStrand; Minhash::Neighbour neighbour;
        queueWrapper.pop(&fragment, &partition, &forwardStrand, &neighbour);
        SamWriter::Alignment retAlignment; int retScore;
        bool happy = alignMinhashNeighbour_long(currentRead, fragment, partition, forwardStrand, &referenceGenomeBrigde,
                                                neighbour, &retAlignment, &score);
        if (retScore > bestScore) {
            bestAlignment = retAlignment;
            bestScore = score;
            foundOne = happy;
        }
    }

    if (bestAlignment.qname.compare("*")==0) {
        makeUnmapped_long(currentRead, bestAlignment);
    }
    return bestAlignment;
}

void align_long(std::string &fastqFile, int &tfBatchSize, TensorflowInference &inferEngine,
                  std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde,
                  SamWriter &samWriter, int nThreads) {
    FastQ fastq(fastqFile);
    int windowLength = referenceGenomeBrigde.getWindowLength();
    int strides = windowLength/2;

    while (fastq.hasNext()) {
        std::vector<LongReadsWrapper > readsVector(tfBatchSize);
        int loadCount = 0;
        int totalFragments = 0;
        for (int i=0; i<tfBatchSize && fastq.hasNext(); i++) {
            LongReadsWrapper readsWrapper;
            readsWrapper.read.reset( new InputRead(fastq.next()) );
            readsWrapper.reverseRead.reset(new std::string(reverseComplement(readsWrapper.read->sequence)));
            int totalKmers = 0;

            // Compute number of numFragments and corresponding starting points
            int end_coord = readsWrapper.read->sequence.length() - windowLength;
            int numFragments = 0;
            int estimateOfFragments = ceil(readsWrapper.read->sequence.length()/strides) + 5;
            int *startCoords = new int[estimateOfFragments];
            if (end_coord >= 0) {
                for (int i = 0, start = 0;; start = std::min(start + strides, end_coord), i++) {
                    numFragments++;
                    startCoords[i] = start;
                    if (start == end_coord)
                        break;
                }
            }
            else {
                startCoords[0]=0;
                numFragments = 1;
            }

            readsWrapper.numFragments.reset(new int(numFragments));
            readsWrapper.positiveFragments.reset(new std::vector<std::string>(numFragments));
            readsWrapper.kmers.reset(new std::vector<std::shared_ptr<Kmer> >(numFragments));
            readsWrapper.revKmers.reset(new std::vector<std::shared_ptr<Kmer> >(numFragments));
            for (int i=0; i<numFragments; i++) {
                int start = startCoords[i];
                std::string positiveWindow = readsWrapper.read->sequence.substr(start, windowLength);
                std::string negativeWindow = reverseComplement(positiveWindow);
                readsWrapper.positiveFragments->at(i) = (positiveWindow);
                readsWrapper.kmers->at(i) = ( std::shared_ptr<Kmer>( encodeWindow(positiveWindow, &totalKmers) ) );
                readsWrapper.revKmers->at(i) = std::shared_ptr<Kmer>( encodeWindow(negativeWindow, &totalKmers) );
            }
            readsWrapper.totalKmers.reset( new int(totalKmers) );
            delete [] startCoords;


            readsVector[i] = readsWrapper;
            loadCount++;
            totalFragments += numFragments;
        }

        prediction_long(inferEngine, readsVector, totalFragments, loadCount, mhIndices.size());

        {
            ThreadPool threadPool(nThreads);
            for (int i = 0; i < loadCount; i++) {
                LongReadsWrapper eachRead = readsVector[i];
                threadPool.enqueue([eachRead, &mhIndices, &referenceGenomeBrigde, &samWriter] {
                    SamWriter::Alignment bestAlignment = alignOne_long(eachRead, mhIndices, referenceGenomeBrigde);
                    samWriter.writeAlignment(bestAlignment);
                });
            }
        }

        LOG(INFO) << loadCount << " finished..";
//        break;
    }
}

#endif