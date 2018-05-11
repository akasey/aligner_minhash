//
// Created by Akash Shrestha on 5/10/18.
//

#include "include/cxxopts.h"
#include "include/common.h"
#include "include/tf_metaparser.h"
#include "include/tensorflow_inference.h"

int main(int argc, const char* argv[]) {
    std::string tfModelDir = "" ;
    std::string mhIndexDir = "" ;
    std::string fastqFile = "" ;
    int nThreads = 4;
    cxxopts::Options options("Aligner", "Does alignment :)");
    options.add_options()
            ("m,tf_dir", "Directory where frozen_graph.pb is located", cxxopts::value<std::string>(tfModelDir), "TF Model dir")
            ("t,threads", "No. of threads", cxxopts::value<int>(nThreads), "Num Threads")
            ("i,minhash_dir", "Directory where all index-xx.mh files are located", cxxopts::value<std::string>(mhIndexDir), "Minhash index dir")
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

    if (tfModelDir.empty() || mhIndexDir.empty() || fastqFile.empty()) {
        std::cerr << options.help() << std::endl;
        exit(-1);
    }

    // @TODO
    TF_MetaParser tfMeta(tfModelDir+"/aligner.meta");
    TensorflowInference inferEngine(tfModelDir+"/frozen_graph.pb", "inputs", "output_logits", tfMeta.getInt("K"), nThreads);
    std::vector<std::string> indices = getFilesInDirectory(mhIndexDir, ".mh");
    std::for_each(indices.begin(), indices.end(), [](std::string i) -> void { std::cout << i << std::endl; });

    return 0;
}
