//
// Created by Akash Shrestha on 5/9/18.
//

#include "../include/common.h"
#include "../include/kmer.h"
#include "../include/tensorflow_inference.h"
#include "../include/ThreadPool.h"

int main() {

    ThreadPool pool(4);
    for (int i=0; i<10; i++) {
        pool.enqueue([i] {
            int k = 7;
            std::string sequence = "TCTGAGCACTTTTGTACTTTGTCATCTGACTAAAAAGGCGTCGAAGCGCCTTTAAAAAATAGTCGAATCAGTAAATTACTGGTATTCGCTAATCGGTACGCAGGCGCAGAACAGGTTACGGTCGCCGTAAACATCATCCAGACGTTCCACTGTCGGCCAGTATTTGTCTGCCACACCTGCCGGGAATACCGCACCTTCAC";

            TensorflowInference inferEngine(
                    "/Users/akash/PycharmProjects/aligner/sample_classification_run/model_dir/tensorboard/frozen_graph.pb",
                    "inputs", "output_logits", 7, 4);

            int totalKmers;
            Kmer *kmers = encodeWindow(sequence, &totalKmers);
            std::set<int> categories = inferEngine.inference(kmers, totalKmers);
            delete[] kmers;

            LOG(INFO) << "From task " << i;
            std::for_each(categories.begin(), categories.end(), [](int i) -> void { std::cout << i << std::endl; });
        });
    }

    return 0;
}