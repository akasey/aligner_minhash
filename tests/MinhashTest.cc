//
// Created by Akash Shrestha on 5/7/18.
//

#include <iostream>
#include <mutex>
#include "../include/cxxopts.h"
#include "../include/ThreadPool.h"
#include "../include/minhash.h"
#include "../include/fastamm.h"
#include "../include/indexer_jobparser.h"

/*
std::map<int, std::string> makeSampleWindows() {
    std::string sequence = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
            "TTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAA"
            "TATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACC"
            "ATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAG"
            "GAGTCCACCCGCCGTATTGCGGCAAGCCGCATTCCGGCTGATCACATGGTGCTGATGGCAGGTTTCACCG"
            "CCGGTAATGAAAAAGGCGAACTGGTGGTGCTTGGACGCAACGGTTCCGACTACTCTGCTGCGGTGCTGGC"
            "TGCCTGTTTACGCGCCGATTGTTGCGAGATTTGGACGGACGTTGACGGGGTCTATACCTGCGACCCGCGT"
            "CAGGTGCCCGATGCGAGGTTGTTGAAGTCGATGTCCTACCAGGAAGCGATGGAGCTTTCCTACTTCGGCG";
    return makeSlidingWindow(sequence, 0, 200, 5);
};


int main(int argc, char *argv[]) {
    Minhash mh;
    std::map<int,std::string> windows = makeSampleWindows();
    for (std::map<int, std::string>::iterator itr=windows.begin(); itr!=windows.end(); itr++) {
        mh.addDocument(itr->first, itr->second);
    }
    FILE *fout = fopen("/Users/shrestha/CSU/aligner_minhash/cmake-build-debug/temp/minhash.index", "wb");
    mh.serialize(fout);
    fclose(fout);

    std::set<Minhash::Neighbour> neighbours = mh.findNeighbours("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCAT");
    for (std::set<Minhash::Neighbour>::iterator itr=neighbours.begin(); itr!=neighbours.end(); itr++) {
        std::cout << "DocID: " << itr->id << ", jaccard: " << itr->jaccard << std::endl;
    }

    Minhash dh;
    FILE *fin = fopen("/Users/shrestha/CSU/aligner_minhash/cmake-build-debug/temp/minhash.index", "rb");
    dh.deserialize(fin);

    dh.compareTest(mh);
}
*/


int createEachMinhashIndex(std::string &baseDirectory, const int segId, std::pair<int, std::string> &segmentPair, int windowLen) {
    int strides = KMER_K -1;
    Minhash mh;
    int startingPoint = segmentPair.first;
    std::unordered_map<int, std::string> slidingWindows = makeSlidingWindow(segmentPair.second, startingPoint, windowLen, strides);
    for (std::unordered_map<int, std::string>::iterator itr=slidingWindows.begin(); itr!=slidingWindows.end(); itr++) {
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
    return 0;
}

void createIndices(std::string &baseDirectory, int nThreads) {
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
        ThreadPool threadPool(nThreads);
        std::vector< std::future<int> > results;
        jobParser.prepareClassificationJob(&fasta);
        IndexerJobParser::Iterator iterator = jobParser.makeJobIterator();
        while(iterator.hasNext()) {
            std::pair<const int, std::pair<int, std::string>> itr = iterator.next();
            int segId = itr.first;
            std::pair<int, std::string> segmentPair = itr.second;
            results.emplace_back(
                    threadPool.enqueue(createEachMinhashIndex, baseDirectory, segId, segmentPair, windowLength)
            );
        }
    }
    catch (ElementNotFoundException &e){
        LOG(ERROR) << e.what();
        exit(-69);
    }
    LOG(INFO) << "Indexing complete...";
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

int main(int argc, const char * argv[]) {
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

    LOG(INFO) << "Reading in minhash indices...";

    std::string mhIndexDir = baseDirectory;
    std::map<std::string, Minhash *> mhIndices;
    std::vector<std::string> indices = getFilesInDirectory(mhIndexDir, ".mh");
    readIndices(indices, nThreads, &mhIndices);

    while (true) {
        std::string mystr;
        std::cout << "Enter sequence: " << std::endl;
        getline(std::cin, mystr);
        std::vector<std::string> splits = split(mystr, "$");
        if (splits.size() == 1) {
            for (int i = 0; i < indices.size(); i++) {
                std::string key = "index-" + std::to_string(i) + ".mh";
                std::set<Minhash::Neighbour> neighbours = mhIndices[key]->findNeighbours(mystr);
                if (neighbours.size() > 0) {
                    std::cout << "Neighbours: " << neighbours.size() << " on segment: " << i << std::endl;
                }

            }
        }
        else if (splits.size() == 2) {
            int splitV = std::stoi(splits[1]);
            if (splitV < mhIndices.size()) {
                std::string key = "index-" + std::to_string(splitV) + ".mh";
                std::set<Minhash::Neighbour> neighbours = mhIndices[key]->findNeighbours(mystr);
                if (neighbours.size() > 0) {
                    std::cout << "Neighbours: " << neighbours.size() << " on segment: " << splitV << std::endl;
                    for (const Minhash::Neighbour &each : neighbours) {
                        std::cout << each.id << " ";
                    }
                    std::cout << std::endl;
                }
            }
        }
        std::cout << "------------------------------" << std::endl;
    }
}
