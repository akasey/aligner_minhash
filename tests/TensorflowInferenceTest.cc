//
// Created by Akash Shrestha on 5/9/18.
//

#include "../include/common.h"
#include "../include/kmer.h"
#include "../include/tf_metaparser.h"
#include "../include/tensorflow_inference.h"
#include "../include/ThreadPool.h"
#include "../include/indexer_jobparser.h"
#include "../include/fastaQ.h"

struct PairedReadsWrapper {
    std::auto_ptr<InputRead> read[2];
    std::auto_ptr<Kmer> kmer[2];
    std::auto_ptr<Kmer> revKmer[2];
    std::auto_ptr<int> totalKmer;
    std::auto_ptr<std::string> reversedRead[2];
    std::auto_ptr<std::set<std::pair<int, bool> > >predictedSegments[2];
};

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


int main() {
    std::string DIRECTORY = "/Users/akash/PycharmProjects/aligner/sample_classification_run/pairedReads/";

    std::map<int, std::pair<int,int>> segmentsRegion;
    FastaMM fasta(DIRECTORY);
    IndexerJobParser refBridge(DIRECTORY);
    refBridge.prepareClassificationJob(&fasta);
    IndexerJobParser::Iterator iterator = refBridge.makeJobIterator();
    while (iterator.hasNext()) {
        std::pair<const int, std::pair<int, std::string>> job = iterator.next();
        int startPoint = job.second.first;
        int endPoint = startPoint + job.second.second.length();
        segmentsRegion[job.first] = std::pair<int,int> (startPoint, endPoint);
    }

    int K = refBridge.getK();
    int tfBatchSize = 10;


    std::string tfModelDir = DIRECTORY + "frozen-graphs";
    TF_MetaParser tfMeta(tfModelDir+"/aligner.meta");
    TensorflowInference inferEngine(
            tfModelDir + "/frozen_graph.pb",
            tfMeta, 4);


    std::string files[] = { DIRECTORY+"/out.bwa.read1.fastq", DIRECTORY+"/out.bwa.read2.fastq" };
    PairedFastQ fastaLoader( files );

    while (fastaLoader.hasNext()) {
        std::vector<PairedReadsWrapper> readsVector(tfBatchSize);
        std::vector< std::pair<Kmer *, int> > pairs(tfBatchSize*2);
        int loadCount = 0;
        for (int i=0; i<tfBatchSize && fastaLoader.hasNext(); i++) {
            PairedReadsWrapper readsWrapper;
            InputRead pairedReads[2];
            fastaLoader.next(pairedReads);
            readsWrapper.read[0] = std::auto_ptr<InputRead>(new InputRead(pairedReads[0]));
            readsWrapper.read[1] = std::auto_ptr<InputRead>(new InputRead(pairedReads[1]));
            int totalKmer = 0;
            Kmer * kmer0 = encodeWindow(readsWrapper.read[0]->sequence, &totalKmer, K);
            readsWrapper.kmer[0] = std::auto_ptr<Kmer>(kmer0);
            Kmer * kmer1 = encodeWindow(readsWrapper.read[1]->sequence, &totalKmer, K);
            readsWrapper.kmer[1] = std::auto_ptr<Kmer>(kmer1);
            readsWrapper.totalKmer = std::auto_ptr<int>(new int(totalKmer));
            readsVector[i] = readsWrapper;
            std::pair<Kmer *, int> onePair(kmer0, totalKmer);
            std::pair<Kmer *, int> twoPair(kmer1, totalKmer);
            pairs[2*i] = onePair;
            pairs[2*i+1] = twoPair;
            loadCount++;
        }

        pairedPrediction(inferEngine, readsVector, pairs, loadCount);

        for (int i=0; i<loadCount; i++) {
            PairedReadsWrapper *currentRead = &readsVector[i];
            std::vector<std::string> trueRange = split(currentRead->read[0]->key, "_");
            for (int j=0; j<2; j++) {
                std::cout << j << "th pair ";
                for (std::set<std::pair<int, bool>>::iterator itr = currentRead->predictedSegments[j]->begin();
                     itr != currentRead->predictedSegments[j]->end(); itr++) {
                    std::cout << itr->first << " ";
                    if (itr->first < segmentsRegion.size()) {
                        auto region = segmentsRegion[itr->first];
                        std::cout << region.first << "-" << region.second << "->" << trueRange[1 + i] << " ";
                        if (std::stoi(trueRange[1 + j]) >= region.first &&
                            std::stoi(trueRange[1 + j]) <= region.second) {
                            std::cout << " hit ";
                        }
                    } else {
                        std::cout << currentRead->read[j]->key << " may be out ";
                    }
                }
                std::cout << "         " ;
            }
            std::cout << std::endl;
        }
    }

    return 0;
}