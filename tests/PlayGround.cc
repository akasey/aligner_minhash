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
#include "edlib.h"

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

void printAlignment(const char* query, const char* target,
                    const unsigned char* alignment, const int alignmentLength,
                    const int position, const EdlibAlignMode modeCode, int eachLineChars) {
    int tIdx = -1;
    int qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != EDLIB_EDOP_INSERT)
                tIdx--;
        }
    }
    for (int start = 0; start < alignmentLength; start += eachLineChars) {
        // target
        printf("T: ");
        int startTIdx;
        for (int j = start; j < start + eachLineChars && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_INSERT)
                printf("-");
            else
                printf("%c", target[++tIdx]);
            if (j == start)
                startTIdx = tIdx;
        }
        printf(" (%d - %d)\n", std::max(startTIdx, 0), tIdx);

        // match / mismatch
        printf("   ");
        for (int j = start; j < start + eachLineChars && j < alignmentLength; j++) {
            printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
        }
        printf("\n");

        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + eachLineChars && j < alignmentLength; j++) {
            if (alignment[j] == EDLIB_EDOP_DELETE)
                printf("-");
            else
                printf("%c", query[++qIdx]);
            if (j == start)
                startQIdx = qIdx;
        }
        printf(" (%d - %d)\n\n", std::max(startQIdx, 0), qIdx);
    }
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


    std::string query = "Akash";
    std::string reference = "BikAsh";
    EdlibAlignTask task = EDLIB_TASK_PATH;
    EdlibAlignMode mode = EDLIB_MODE_HW;

    char lowercase[] = {'a','c','g','t'};
    char uppercase[] = {'A', 'C', 'G', 'T'};
    EdlibEqualityPair pairs[4];
    for (int i=0; i<4; i++) {
        EdlibEqualityPair pair;
        pair.first = lowercase[i];
        pair.second = uppercase[i];
        pairs[i] = pair;
    }

    EdlibAlignResult result = edlibAlign(query.c_str(), query.length(), reference.c_str(), reference.length(),
               edlibNewAlignConfig(-1, mode, task, pairs, 4));
    if (result.status == EDLIB_STATUS_OK) {
        char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_EXTENDED);
        std::string cigarString(cigar);
        std::cout << cigarString << std::endl;
        std::cout << "EditDistance: " << result.editDistance << std::endl;
        std::cout << "Start: " << *(result.startLocations) << std::endl;
        std::cout << "End: " << *(result.endLocations) << std::endl;
        printAlignment(query.c_str(), reference.c_str(), result.alignment, result.alignmentLength, *(result.endLocations), mode, 200);
    }
    edlibFreeAlignResult(result);
    return 0;

}