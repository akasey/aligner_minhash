//
// Created by Akash Shrestha on 5/10/18.
//

#include "include/sam_writer.h"
#include "include/cxxopts.h"
#include "include/common.h"
#include "include/minhash.h"
#include "include/tf_metaparser.h"
#include "include/tensorflow_inference.h"
#include "include/ThreadPool.h"
#include "include/sam_writer.h"
#include "include/indexer_jobparser.h"
#include "include/priorityqueue_wrapper.h"

void align_single(std::string &fastqFile, int &tfBatchSize, TensorflowInference &inferEngine,
                  std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde,
                  SamWriter &samWriter);
void align_paired(std::string &fastqFiles, int &tfBatchSize, TensorflowInference &inferEngine,
                  std::map<std::string, Minhash *> &mhIndices, IndexerJobParser &referenceGenomeBrigde,
                  SamWriter &samWriter);

std::map<std::string, std::string> processArguments(int argc, const char *argv[]) {
    std::string tfModelDir = "" ;
    std::string mhIndexDir = "" ;
    std::string referenceGenomeDir = "";
    std::string fastqFile = "" ;
    std::string samFile = "";
    std::string mode = "single";
    int nThreads = 4;
    int tfBatchSize = 5;

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
            ("mode", "Single or Paired alignment", cxxopts::value<std::string>(mode), "[single|paired]")
            ("o,output", "Output SamFile", cxxopts::value<std::string>(samFile), "Sam file");


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

    if (mode.compare("single") != 0 && mode.compare("paired") != 0) {
        std::cerr << "mode: [single|paired]\n";
        std::cerr << options.help() << std::endl;
        exit(-1);
    }

    if (samFile.empty()) {
        samFile = referenceGenomeDir + "/alignment.sam";
    }

    typedef std::pair<std::string,std::string> __;
    std::map<std::string, std::string> map = {
            __("tensorflow", tfModelDir),
            __("minhash", mhIndexDir),
            __("reference", referenceGenomeDir),
            __("fastq", fastqFile),
            __("output", samFile),
            __("batchsize", std::to_string(tfBatchSize)),
            __("threads", std::to_string(nThreads)),
            __("mode", mode),
            __("command", wholeCommand)

    };
    return map;
}

#if THREADS_ENABLE
void readIndices(std::vector<std::string> mhIndexLocations, int nThreads, std::map<std::string, Minhash *> *mhIndices) {
    std::map<std::string, std::future<Minhash *> > mhThreadResults;
    {
        ThreadPool threadPool(nThreads);
        for (int i = 0; i < mhIndexLocations.size(); i++) {
            std::string filename = mhIndexLocations[i];
            std::string basefname = basename(filename);
            mhThreadResults[basefname] = threadPool.enqueue([i, filename] {
                Minhash *mh = new Minhash();
                mh->deserialize(filename);
                return mh;
            });
        }
    } // end of parallelization

    for (std::map<std::string, std::future<Minhash *> >::iterator itr=mhThreadResults.begin(); itr!=mhThreadResults.end(); itr++) {
        (*mhIndices)[itr->first] = itr->second.get();
    }
}
#else
void readIndices(std::vector<std::string> mhIndexLocations, int nThreads, std::map<std::string, Minhash *> *mhIndices) {
    for (int i = 0; i < mhIndexLocations.size(); i++) {
        std::string filename = mhIndexLocations[i];
        std::string basefname = basename(filename);
        Minhash *mh = new Minhash();
        mh->deserialize(filename);
        (*mhIndices)[basefname] = mh;
    }
}
#endif

int main(int argc, const char* argv[]) {
    std::map<std::string, std::string> arguments = processArguments(argc, argv);
    int nThreads = stoi(arguments["threads"]);
    int tfBatchSize = stoi(arguments["batchsize"]);

    // Actual main
    IndexerJobParser referenceGenomeBrigde(arguments["reference"]);
    SamWriter samWriter(SamWriter::Mode::SINGLE, arguments["output"]);
    { // scope to free some memory
        FastaMM fastamm(arguments["reference"]);
        referenceGenomeBrigde.prepareClassificationJob(&fastamm);
        samWriter.writeHeaders(arguments["command"], fastamm);
    }



    LOG(INFO) << "Reading tensorflow graph...";
    TF_MetaParser tfMeta(arguments["tensorflow"]+"/aligner.meta");
    TensorflowInference inferEngine(arguments["tensorflow"]+"/frozen_graph.pb", tfMeta, nThreads);
    std::vector<std::string> indices = getFilesInDirectory(arguments["minhash"], ".mh");
    int K = tfMeta.getInt("K");
    int numClasses = tfMeta.getInt("output_shape");
    if (numClasses-1 != indices.size()) {
        LOG(ERROR) << "Minhash indices doesn't equals to neural network classes";
        exit(-1);
    }


    LOG(INFO) << "Reading in minhash indices...";
    std::map<std::string, Minhash *> mhIndices;
    readIndices(indices, nThreads, &mhIndices);

    LOG(INFO) << "KMER_K " << KMER_K;

/*    for (int i =0; i<indices.size(); i++) {
        std::string base = basename(indices[i]);
        std::cout << mhIndices[base]-> filename  << " " << base << std::endl;
    }*/

    if (arguments["mode"].compare("single") == 0) {
        align_single(arguments["fastq"], tfBatchSize, inferEngine, mhIndices, referenceGenomeBrigde, samWriter);
    }
    else if (arguments["mode"].compare("paired") == 0) {
        align_paired(arguments["fastq"], tfBatchSize, inferEngine, mhIndices, referenceGenomeBrigde, samWriter);
    }

    for (std::map<std::string, Minhash *>::iterator itr = mhIndices.begin(); itr != mhIndices.end(); itr++) {
        delete itr->second;
    }
    return 0;
}
