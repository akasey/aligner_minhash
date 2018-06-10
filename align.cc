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

int main(int argc, const char* argv[]) {
    std::string tfModelDir = "" ;
    std::string mhIndexDir = "" ;
    std::string referenceGenomeDir = "";
    std::string fastqFile = "" ;
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
            ("f,fastq", "Input FastQ file for aligning", cxxopts::value<std::string>(fastqFile), "FastQ file");


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

    // Actual main
    FastaMM fastamm(referenceGenomeDir);
    SamWriter samWriter(SamWriter::Mode::SINGLE, "x.sam");
    samWriter.writeHeaders(wholeCommand, fastamm);


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
        std::map<std::string, std::pair<Kmer *, int> > readsMap;
        std::vector< std::pair<Kmer *, int> > pairs(tfBatchSize);
        for (int i=0; i<tfBatchSize && fastq.hasNext(); i++) {
            InputRead read = fastq.next();
            int totalKmers;
            Kmer *kmers = encodeWindow(read.sequence, &totalKmers);
            std::pair<Kmer *, int> onePair(kmers, totalKmers);
            readsMap[read.key] = onePair;
            pairs[i] = onePair;
        }

        Tensor tensor = inferEngine.makeTensor(pairs);
        std::vector< std::set<int> > predictions = inferEngine.inference(tensor);

        int i=0;
        for (std::map<std::string, std::pair<Kmer *, int> >::iterator itr=readsMap.begin(); itr != readsMap.end(); itr++) {
            std::string readKey = itr->first;
            std::set<int> predicted = predictions[i];
            for (std::set<int>::iterator itrr = predicted.begin(); itrr != predicted.end(); itrr++) {
                std::string key = "index-" + std::to_string(*itrr) + ".mh";
                if (*itrr < mhIndices.size()) {
//                    std::set<Minhash::Neighbour> neighbours = mhIndices[key]->findNeighbours(itr->second.first, itr->second.second);
                    std::set<Minhash::Neighbour> neighbours = mhIndices[key]->findNeighbours(itr->second.first, itr->second.second);
                    std::cout << "Predictions for " << readKey << " key: " << key << " size: " << neighbours.size() << std::endl;
                    std::for_each(neighbours.begin(), neighbours.end(), [](Minhash::Neighbour a) {
                        std::cout << a.id << "(" << a.jaccard << "), ";
                    });
                    std::cout << std::endl;
                }
                else
                    LOG(ERROR) << readKey << " predicted as " << std::to_string(*itrr);
            }
            i++;

            delete [] readsMap[readKey].first;
        }
    }


    for (std::map<std::string, Minhash *>::iterator itr = mhIndices.begin(); itr != mhIndices.end(); itr++) {
        delete itr->second;
    }
    return 0;
}
