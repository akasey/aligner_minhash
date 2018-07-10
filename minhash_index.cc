//
// Created by Akash Shrestha on 5/8/18.
//

#include "include/cxxopts.h"
#include "include/common.h"
#include "include/fastamm.h"
#include "include/indexer_jobparser.h"
#include "include/minhash.h"
#include "include/ThreadPool.h"

void createEachMinhashIndex(std::string baseDirectory, const int segId, std::pair<int, std::string> segmentPair, int windowLen) {
    int strides = KMER_K -1;
    int startingPoint = segmentPair.first;
    std::map<int, std::string> slidingWindows = makeSlidingWindow(segmentPair.second, startingPoint, windowLen, strides);
    Minhash mh;
    for (std::map<int, std::string>::iterator itr=slidingWindows.begin(); itr!=slidingWindows.end(); itr++) {
        mh.addDocument(itr->first, itr->second);
    }
    std::string filename = baseDirectory +"/indices/index-" + std::to_string(segId) + ".mh";
    FILE *stream = fopen(filename.c_str(), "wb");
    if (!stream) {
        LOG(ERROR) << "Can't create stream.. LET ME DIE!" << filename;
        LOG(INFO) << "Trying to create directory for indices...";
        mkdir(baseDirectory+"/indices");
        stream = fopen(filename.c_str(), "wb");
        if (!stream) {
            LOG(ERROR) << "Really really can't.. PLEASE LET ME DIE!!!" << filename;
            exit(-2);
        }
    }
    mh.serialize(stream);
    fclose(stream);
    LOG(INFO) << segId << " indexing complete..";
}

void mainStuffs(std::string &baseDirectory, int nThreads) {
    LOG(INFO) << "Basedir: " << baseDirectory;
    LOG(INFO) << "Reading " << baseDirectory + "/sequence.fasta...";
    FastaMM fasta(baseDirectory);
    LOG(INFO) << "Reading " << baseDirectory + "/classify_detail.log...";
    IndexerJobParser jobParser(baseDirectory);
    int numClasses = jobParser.getNumClasses();
    NULL_CHECK(numClasses, "Number of classes not in " + baseDirectory);
    LOG(INFO) << "Need to index " << numClasses << " indices....";

    int windowLength = jobParser.getWindowLength();
    NULL_CHECK(windowLength, "WindowLength not in " + baseDirectory);
    LOG(INFO) << "WindowLength " << windowLength << "....";

    try {
#if THREADS_ENABLE
        ThreadPool threadPool(nThreads);
#endif
        jobParser.prepareClassificationJob(&fasta);
        IndexerJobParser::Iterator iterator = jobParser.makeJobIterator();
        while(iterator.hasNext()) {
            std::pair<const int, std::pair<int, std::string>> itr = iterator.next();
            int segId = itr.first;
            std::pair<int, std::string> segmentPair = itr.second;
#if THREADS_ENABLE
            threadPool.enqueue([baseDirectory, segId, segmentPair, windowLength] {
                createEachMinhashIndex(baseDirectory, segId, segmentPair, windowLength);
            });
#else
            createEachMinhashIndex(baseDirectory, segId, segmentPair, windowLength);
#endif
        }
    }
    catch (ElementNotFoundException &e){
        LOG(ERROR) << e.what();
        exit(-69);
    }
    LOG(INFO) << "Indexing complete...";
}


int main(int argc, const char* argv[]) {
    std::string baseDirectory = "" ;
    int nThreads = 5;
    cxxopts::Options options("Minhash Indexer", "Create the minhash index for after classification.");
    options.add_options()
            ("d,dir", "Directory where sequence.fasta and classify_detail.log are placed.", cxxopts::value<std::string>(baseDirectory)->default_value("./"), "Base directory")
            ("t,threads", "No. of threads", cxxopts::value<int>(nThreads), "Num Threads");
    if (argc <= 1){
        std::cerr << options.help() << std::endl;
        exit(-1);
    }

    auto result = options.parse(argc, argv);
    mainStuffs(baseDirectory, nThreads);
    return 0;
}