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

struct PairedReadsWrapper {
    InputRead *read[2];
    Kmer *kmer[2];
    int *totalKmer;
    std::string *reversedRead[2];
    std::set<std::pair<int, bool>> *predictedSegments[2];

    void clear() {
        if (totalKmer) delete totalKmer;
        for (int i=0; i<2; i++) {
            if (read[i]) delete[] read[i];
            if (kmer[i]) delete[] kmer[i];
            if (reversedRead[i]) delete reversedRead[i];
            if (predictedSegments[i]) delete predictedSegments[i];
        }
    }
};


void pairedPrediction(TensorflowInference &inferEngine, std::vector<PairedReadsWrapper> &readsVector, std::vector< std::pair<Kmer*,int> > &pairs, int loadCount) {
    Tensor tensor = inferEngine.makeTensor(pairs); // pairs.size() = 2*loadCount
    std::vector<std::set<std::pair<int, bool> > > predictions = inferEngine.inference(tensor);
    for (int i=0; i<loadCount; i++) {
        readsVector[i].predictedSegments[0] = new std::set<std::pair<int, bool> >();
        readsVector[i].predictedSegments[1] = new std::set<std::pair<int, bool> >();
        std::set<std::pair<int, bool> > inner = predictions[i];
        for (int j=0; j<predictions.size(); j++) {
            std::pair<int, bool> lacasito = predictions[j];
            int key = j%2;
            readsVector[i].predictedSegments[key]->insert(lacasito);
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
    }
}