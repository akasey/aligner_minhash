//
// Created by Akash Shrestha on 5/9/18.
//

#ifndef ALIGNER_MINHASH_TENSORFLOW_INFERENCE_H
#define ALIGNER_MINHASH_TENSORFLOW_INFERENCE_H

#include "common.h"
#include "kmer.h"
#include "tf_metaparser.h"
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
    TF_MetaParser tfMeta;
    string input_layer;
    string output_layer;
    int numCore;
    int K;
    int input_shape, output_shape;
    std::unique_ptr<tensorflow::Session> session;

    Status LoadGraph(const string& graph_file_name, std::unique_ptr<tensorflow::Session>* session);
//    Status GetIOShapes(tensorflow::GraphDef &graphdef);
    Status GetNonZeroLabels(const std::vector<Tensor>& outputs, Tensor* coords);

public:
    TensorflowInference(std::string inputGraph, TF_MetaParser &tf_meta, int numThreads);
    ~TensorflowInference();
//    TensorflowInference(const TensorflowInference &p2); // copy constructor

    std::vector<std::set<int> > inference(Tensor);
    Tensor makeTensor(std::vector< std::pair<Kmer *, int> > pairs);
    int getOutputShape();
    int getInputShape();
    int getK();
};


#endif //ALIGNER_MINHASH_TENSORFLOW_INFERENCE_H
