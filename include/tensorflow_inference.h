//
// Created by Akash Shrestha on 5/9/18.
//

#ifndef ALIGNER_MINHASH_TENSORFLOW_INFERENCE_H
#define ALIGNER_MINHASH_TENSORFLOW_INFERENCE_H

#include "common.h"
#include "kmer.h"
#include "tensorflow/core/platform/env.h"
#include "tensorflow/core/public/session.h"
#include "tensorflow/core/util/command_line_flags.h"
#include "tensorflow/core/platform/init_main.h"
#include "tensorflow/core/lib/io/path.h"
#include "tensorflow/cc/framework/scope.h"
#include "tensorflow/cc/ops/standard_ops.h"

using tensorflow::Flag;
using tensorflow::Tensor;
using tensorflow::Status;
using tensorflow::string;
using tensorflow::int64;

class TensorflowInference {
private:
    std::string inputGraph;
    string input_layer;
    string output_layer;
    int numCore;
    int K;
    int input_shape, output_shape;
    std::unique_ptr<tensorflow::Session> session;

    Status LoadGraph(const string& graph_file_name, std::unique_ptr<tensorflow::Session>* session);
    Status GetIOShapes(tensorflow::GraphDef &graphdef);
    Status GetNonZeroLabels(const std::vector<Tensor>& outputs, Tensor* coords);
    Tensor makeTensor(Kmer *kmers, int totalKmers);

public:
    TensorflowInference(std::string inputGraph, std::string ip, std::string op, int K, int numThreads);
    ~TensorflowInference();
//    TensorflowInference(const TensorflowInference &p2); // copy constructor

    std::set<int> inference(Kmer *kmers, int totalKmers);
    int getOutputShape();
    int getInputShape();
};


#endif //ALIGNER_MINHASH_TENSORFLOW_INFERENCE_H
