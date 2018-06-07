//
// Created by Akash Shrestha on 5/9/18.
//

#include "../include/common.h"
#include "../include/kmer.h"
#include "../include/tf_metaparser.h"
#include "../include/tensorflow_inference.h"
#include "../include/ThreadPool.h"

int main() {

    ThreadPool pool(4);
    for (int i=0; i<1; i++) {
        pool.enqueue([i] {
            int k = 7;
            std::string sequence[] = //"AGGAATTTTTGAAAAAGTACACAGATAAAATCTGATGTCTTTTACGCATAGGTTTGAATCTACCAGCCTGCGTGAATACGATATTCGCGGAATTGTTGGTAAAACCCTTTCAACCGATGATGCGTTTGCCATTGGCCGCACCTTTGGCAGCCTTGTGCGCCGCAACGGCGGTAAAACGCTGGATGTGGGCTATGATGGCC";
                        {"GTTCGTGTTGGCATGTCTGTTGAAGTGAATATTGATACGGCGTCAAAGCCAGAAGGCCCGCATGCCAATAATCCACGTTATATCTGACAGTAATGTCATGAAGCGTCTTCTTCTTTCCGGCTTATCTGTGCTGGTGCTTTCTGGCTGCACAGTCGGGCCAAGCTATCATAAGCCGCAGGTTAAAGCGCCAGCAGATTGGC",
                         "ATGCCAGCAATACGGGCCGTTGCCAGAACATCACCTTTTGAGGCTGCGCCTGAAAGGATCATTTCCAGTGTTTGGGGCTGCATAATAATAGCCCCGGTGGCAACAGCCTTGCGGTGTGTGGCGGCTCTGGCGGAAACATACACCATATGCGCATGCCCCGCCGCATCCAGATGTGTGAGTGCTGGAGAGGGCACGTTAGG" };

            std::string tfModelDir = "/Users/akash/PycharmProjects/aligner/sample_classification_run/model_dir/frozen-graphs";
            TF_MetaParser tfMeta(tfModelDir+"/aligner.meta");
            TensorflowInference inferEngine(
                    tfModelDir + "/frozen_graph.pb",
                    tfMeta, 4);

            std::vector< std::pair<Kmer *, int> > kmerVector(2);
            for (int j=0;j < 2; j++) {
                int totalKmers;
                Kmer *kmers = encodeWindow(sequence[j], &totalKmers);
                kmerVector[i] = std::pair<Kmer *, int>(kmers, totalKmers);
            }
            Tensor tensor = inferEngine.makeTensor(kmerVector);
            std::vector<std::set<int> > categories = inferEngine.inference(tensor);

            LOG(INFO) << "From task " << i;
            for (int j=0; j < categories.size(); j++) {
                std::set<int> predictions = categories[j];
                for (std::set<int>::iterator itr = predictions.begin(); itr != predictions.end(); itr++) {
                    std::cout << j << " " << *itr << std::endl;
                }
            }

            for (int j=0; j<kmerVector.size(); j++) {
                delete[] kmerVector[j].first;
            }
        });
    }

    return 0;
}