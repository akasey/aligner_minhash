//
// Created by Akash Shrestha on 6/11/18.
//

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstdlib>
#include <algorithm>
#include "tf_metaparser.h"
#include "priorityqueue_wrapper.h"
#include "tensorflow_inference.h"

std::set<Minhash::Neighbour> makeRandom() {
    std::set<Minhash::Neighbour> neigh;
    int max = rand() % 10 + 1;
    for (int i=0; i< max; i++) {
        Minhash::Neighbour n;
        n.id = i+1;
        n.jaccard = rand() % 10;
        neigh.insert(n);
        std::cout << "NewNeighbour: (" << n.id << "," << n.jaccard << ")\n";
    }
    return neigh;
}

//	Align a pair of genome sequences.
int main (int argc, char * const argv[]) {
//    PriorityQueueWrapper wrapper;
//    srand(123);
//    int max = 5;
//    for (int i=0; i<max; i++) {
//        std::set<Minhash::Neighbour> queue1 = makeRandom();
//        wrapper.addQueue(i, i%2, &queue1);
//    }
//
//
//    while(wrapper.hasNext()) {
//        Minhash::Neighbour neighbour;
//        int partition;
//        bool forwardStrand;
//        wrapper.pop(&partition, &forwardStrand, &neighbour);
//        std::cout << "partition: " << partition << " forwardStrand: " << forwardStrand << " neighbour:(" << neighbour.id << " ," << neighbour.jaccard  << ")" << std::endl;
//    }
/*
    std::string directory = "/Users/akash/ClionProjects/aligner_minhash/cmake-build-debug/ecoli/model_dir/frozen-graphs/";
    TF_MetaParser meta (directory + "aligner.meta");
    TensorflowInference inference(directory + "frozen_graph.pb", meta, 1);
    int K = meta.getInt("K");

    std::string sequence = "TGCTGCCACGCGTGAGCGGTCGTAATCAGCACCGCATCAGCAAGTGTATCTGCCGTGCACTGCAACAACGCTGCTTCGGCCTGGTAATGGCCCGCCGCCTTCCAGCGTTCGACCCAGGCGTTAGGGTCAATGCGGGTCGCTTCACTTACGCCAATGTCGTTATCCAGCGGTGCACGGGTGAACTGATCGCGCAGCGGCGT";
    int totalKmers = 0;
    std::shared_ptr<Kmer> kmers = std::shared_ptr<Kmer>(encodeWindow(sequence, &totalKmers, K));

    std::vector<std::pair<std::shared_ptr<Kmer>, int> > pairs;
    pairs.push_back({kmers, totalKmers});
    Tensor tensor = inference.makeTensor(pairs);
    std::vector<std::set<std::pair<int, bool> > > predictions = inference.inference(tensor);

    std::vector<Kmer> allKmers;
    for (int i=0; i< totalKmers; i++) {
        allKmers.push_back(kmers.get()[i]);
    }

    std::sort(allKmers.begin(), allKmers.end());
    std::cout << "Total kmers " << totalKmers << std::endl;
    for (int i=0; i<totalKmers; i++) {
        std::cout << allKmers[i] << "\n";
    }
    std::cout << std::endl;

    for (int i=0; i < predictions.size(); i++) {
        std::cout << "Predictions: " << std::endl;
        for (const auto inner : predictions[i]) {
            std::cout << inner.first << " " << inner.second << std::endl;
        }
    }
    */
    srand(10);
    auto MAX_UINT = std::numeric_limits<Kmer>::max();
    for (int i=0; i<30; i++) {
        std::cout << rand() % MAX_UINT << " ";
    }
    std::cout << std::endl;

    uint32_t seed = 10;
    for (int i=0; i<30; i++) {
        std::cout << rand_r(&seed) % MAX_UINT << " ";
    }
    std::cout << std::endl;


}