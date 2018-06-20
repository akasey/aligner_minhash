//
// Created by Akash Shrestha on 5/10/18.
//

#include <sam_writer.h>
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


void processArguments(int argc, const char *argv[]) {

}

std::map<std::string, Minhash *> readIndices(std::vector<std::string> mhIndexLocations, int nThreads) {
    std::map<std::string, Minhash *> mhIndices;
    ThreadPool threadPool(nThreads);
    for (int i = 0; i < mhIndexLocations.size(); i++) {
        std::string filename = mhIndexLocations[i];
        std::string basefname = basename(filename);
        mhIndices[basefname] = new Minhash();
        threadPool.enqueue([i, filename, &mhIndices, basefname] {
            mhIndices[basefname]->deserialize(filename);
        });
    }
    return mhIndices;
}

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

inline bool tryFirstOutofGiven(std::string &queryString, std::set<Minhash::Neighbour> &neighbours,
                               std::pair<int, std::string> *referenceSegment, SamWriter &samWriter,SamWriter::Alignment *retAlignment, int *score) {
    if (neighbours.size() > 0) {
        Minhash::Neighbour first = *(neighbours.begin());
        neighbours.erase(neighbours.begin());
        int start = fmax((int)(first.id) - (int)(referenceSegment->first) - 10, 0);
        int length = fmin(referenceSegment->second.length(), start+220) - start;
        std::string partOfReference = referenceSegment->second.substr(start, length);
        *score = samWriter.alignment(partOfReference, queryString, retAlignment);
        if ( queryString.length()*0.8 <= *score ) { // atleast 80% matches
            return true;
        }
    }
    return false; // not happy with this alignment
}

int main(int argc, const char* argv[]) {
    std::string tfModelDir = "" ;
    std::string mhIndexDir = "" ;
    std::string referenceGenomeDir = "";
    std::string fastqFile = "" ;
    std::string samFile = "";
    int nThreads = 4;
    int tfBatchSize = 1;

    std::string wholeCommand = "";
    for (int i=0; i<argc; i++) {
        wholeCommand += std::string(argv[i]) + " ";
    }

    cxxopts::Options options("Aligner", "Does alignment :)");
    options.add_options()
            ("m,tf_dir", "Directory where frozen_graph.pb is located", cxxopts::value<std::string>(tfModelDir), "TF Model dir")
            ("t,threads", "No. of threads", cxxopts::value<int>(nThreads), "Num Threads")
            ("i,minhash_dir", "Directory where all index-xx.mh files are located", cxxopts::value<std::string>(mhIndexDir), "Minhash index dir")
            ("r,genome_dir", "Directory where sequence.fasta, classify_detail.log are located", cxxopts::value<std::string>(referenceGenomeDir), "Reference genome dir")
            ("f,fastq", "Input FastQ file for aligning", cxxopts::value<std::string>(fastqFile), "FastQ file")
            ("o,output", "Output SamFile", cxxopts::value<std::string>(fastqFile), "Sam file");


    std::vector<cxxopts::KeyValue> arguments;
    try {
        auto result = options.parse(argc, argv);
        arguments = result.arguments();
    }
    catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << options.help() << std::endl;
        exit(-1);
    }

    if (tfModelDir.empty() || mhIndexDir.empty() || fastqFile.empty() || referenceGenomeDir.empty()) {
        std::cerr << options.help() << std::endl;
        exit(-1);
    }
    if (samFile.empty()) {
        samFile = referenceGenomeDir + "/alignment.sam";
    }

    // Actual main
    IndexerJobParser referenceGenomeBrigde(referenceGenomeDir);
    SamWriter samWriter(SamWriter::Mode::SINGLE, samFile);
    { // scope to free some memory
        FastaMM fastamm(referenceGenomeDir);
        referenceGenomeBrigde.prepareClassificationJob(&fastamm);
        samWriter.writeHeaders(wholeCommand, fastamm);
    }



    LOG(INFO) << "Reading tensorflow graph...";
    TF_MetaParser tfMeta(tfModelDir+"/aligner.meta");
    TensorflowInference inferEngine(tfModelDir+"/frozen_graph.pb", tfMeta, nThreads);
    std::vector<std::string> indices = getFilesInDirectory(mhIndexDir, ".mh");
    int K = tfMeta.getInt("K");
    int numClasses = tfMeta.getInt("output_shape");
    if (numClasses-1 != indices.size()) {
        LOG(ERROR) << "Minhash indices doesn't equals to neural network classes";
        exit(-1);
    }


    LOG(INFO) << "Reading in minhash indices...";
    std::map<std::string, Minhash *> mhIndices = readIndices(indices, nThreads);

/*    for (int i =0; i<indices.size(); i++) {
        std::string base = basename(indices[i]);
        std::cout << mhIndices[base]-> filename  << " " << base << std::endl;
    }*/


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
                    std::pair<int, std::string> *referenceSegment = referenceGenomeBrigde.getSegmentForID(pair.first);
                    NULL_CHECK(referenceSegment, "Reference for segment " + std::to_string(pair.first) + " is NULL");

                    std::set<Minhash::Neighbour> posNeighboursCurrentPred = mhIndices[key]->findNeighbours(currentRead->kmer, *(currentRead->totalKmers));

                    SamWriter::Alignment alignment;
                    bool happy = tryFirstOutofGiven(currentRead->read->sequence, posNeighboursCurrentPred, referenceSegment, samWriter, &alignment, &score);
                    if (!happy) {
                        queueWrapper.addQueue(pair.first, true, &posNeighboursCurrentPred);
                    }
                    else {
                        if (score > bestScore) {
                            bestAlignment = alignment;
                            bestScore = score;
                        }
                        continue;
                    }


                    // Also try first of reverse strand
                    std::string reverse = reverseComplement(currentRead->read->sequence);
                    currentRead->reverseRead = new std::string(reverse);
                    int totalKmers = 0;
                    currentRead->revKmer = encodeWindow(reverse, &totalKmers);

                    SamWriter::Alignment negAlignment;
                    std::set<Minhash::Neighbour> negNeighboursCurrentPred = mhIndices[key]->findNeighbours(currentRead->revKmer, *(currentRead->totalKmers));
                    happy = tryFirstOutofGiven(reverse, negNeighboursCurrentPred, referenceSegment, samWriter, &negAlignment, &score);
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
                else
                    LOG(ERROR) << currentRead->read->key << " predicted as " << std::to_string(pair.first);
            }

            while(queueWrapper.hasNext()) {
                int partition; bool forwardStrand; Minhash::Neighbour neighbour;
                queueWrapper.pop(&partition, &forwardStrand, &neighbour);

            }
        }
    }


    for (std::map<std::string, Minhash *>::iterator itr = mhIndices.begin(); itr != mhIndices.end(); itr++) {
        delete itr->second;
    }
    return 0;
}
