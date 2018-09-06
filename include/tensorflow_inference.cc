//
// Created by Akash Shrestha on 5/9/18.
//

#include "tensorflow_inference.h"


TensorflowInference::TensorflowInference(std::string graph_pb, TF_MetaParser &tf_meta, int numThreads=1) : tfMeta(tf_meta) {
    inputGraph = graph_pb;
    input_layer = tfMeta["input_layer"];
    output_layer = tfMeta["output_layer"];
    input_shape = tfMeta.getInt("input_shape");
    output_shape = tfMeta.getInt("output_shape");
    K = tfMeta.getInt("K");
    numCore = numThreads;

    // First we load and initialize the model.
    Status load_graph_status = LoadGraph(inputGraph, &session);
    if (!load_graph_status.ok()) {
        LOG(ERROR) << load_graph_status.ToString();
        throw TensorflowInferenceException("Couldn't load graph");
    }
}

TensorflowInference::~TensorflowInference() {
}

// Reads a model graph definition from disk, and creates a session object you
// can use to run it.
Status TensorflowInference::LoadGraph(const string& graph_file_name, std::unique_ptr<tensorflow::Session>* session) {
    tensorflow::GraphDef graph_def;
    Status load_graph_status =
            ReadBinaryProto(tensorflow::Env::Default(), graph_file_name, &graph_def);
    if (!load_graph_status.ok()) {
        return tensorflow::errors::NotFound("Failed to load compute graph at '",
                                            graph_file_name, "'");
    }

    // initialize the number of worker threads
    tensorflow::SessionOptions options;
    tensorflow::ConfigProto & config = options.config;
    if (numCore > 1) {
        config.set_inter_op_parallelism_threads(numCore);
        config.set_intra_op_parallelism_threads(numCore);
        config.set_use_per_session_threads(false);
    }
    session->reset(tensorflow::NewSession(options));
    Status session_create_status = (*session)->Create(graph_def);
    if (!session_create_status.ok()) {
        return session_create_status;
    }
//    GetIOShapes(graph_def);
    return Status::OK();
}

/*Status TensorflowInference::GetIOShapes(tensorflow::GraphDef &graphdef) {
    // Prepare dummy input
    Kmer *kmers; int totalKmer = 0;
    input_shape = std::pow(4, K);
    Tensor tensor = makeTensor(kmers, totalKmer);
    std::vector<Tensor> outputs;
    Status run_status = session->Run({{input_layer, tensor}},
                                     {output_layer}, {}, &outputs);
    TF_RETURN_IF_ERROR(run_status);

    output_shape = outputs[0].shape().dim_size(1);
    return Status::OK();
}*/

Status TensorflowInference::GetNonZeroLabels(const std::vector<Tensor>& outputs, Tensor* coords) {
    auto root = tensorflow::Scope::NewRootScope();
    using namespace ::tensorflow::ops;  // NOLINT(build/namespaces)

    string output_name = "where_1";
    tensorflow::Output rounded = Round(root.WithOpName("rounder"), outputs[0]).y;
    Where(root.WithOpName(output_name), rounded);
    // This runs the GraphDef network definition that we've just constructed, and
    // returns the results in the output tensors.
    tensorflow::GraphDef graph;
    TF_RETURN_IF_ERROR(root.ToGraphDef(&graph));

    std::unique_ptr<tensorflow::Session> session(
            tensorflow::NewSession(tensorflow::SessionOptions()));
    TF_RETURN_IF_ERROR(session->Create(graph));
    // The TopK node returns two outputs, the scores and their original indices,
    // so we have to append :0 and :1 to specify them both.
    std::vector<Tensor> out_tensors;
    TF_RETURN_IF_ERROR(session->Run({}, {output_name},
                                    {}, &out_tensors));
    *coords = out_tensors[0];
    return Status::OK();
}

Tensor TensorflowInference::makeTensor(std::vector< std::pair<std::shared_ptr<Kmer>, int> > pairs) {
    int totalRows = pairs.size();
    Tensor tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape( { totalRows,input_shape} ));
    for (int i=0; i < totalRows; i++) {
        std::pair<std::shared_ptr<Kmer>, int> kmerRows = pairs[i];
        auto vector = tensor.matrix<float>();
        for (int j = 0; j < input_shape; j++) {
            vector(i, j) = 0;
//            vector(i+1, j) = 0; // reverse strand
//            if ( j==input_shape-1 )
//                vector(i+1, j) = 1; // reverse strand

        }
        for (int j = 0; j < kmerRows.second; j++) {
            Kmer enumeration = kmerRows.first.get()[j];
            vector(i, enumeration) = 1;
//            vector(i+1, enumeration) = 1;
        }
    }
    return tensor;
}

std::vector<std::set<std::pair<int, bool> > > TensorflowInference::inference(Tensor tensor) {
    int numRows = tensor.dim_size(0);
    std::vector<Tensor> outputs;
    Status run_status = session->Run({{input_layer, tensor}},
                                     {output_layer}, {}, &outputs);
    if (!run_status.ok()) {
        LOG(ERROR) << "Running model failed: " << run_status;
        throw TensorflowInferenceException("Session::Run failed..");
    }

    Tensor coord;
    run_status = GetNonZeroLabels(outputs, &coord);
    if (!run_status.ok()) {
        LOG(ERROR) << "GetNonZeroLabels failed " << run_status;
        throw TensorflowInferenceException("GetNonZeroLabels failed..");
    }

    std::vector<std::set<std::pair<int, bool> > > toReturn(numRows);
//    std::set<int> predictions;
    for (int i=0; i<coord.NumElements(); i+=2) {
        int rowNum = coord.matrix<int64>()(i);
        int category = coord.matrix<int64>()(i+1);
//        if (predictions.find(category) == predictions.end()) {
            // As even rowNums are predictions for positive strand. We discard prediction of negative if positive also belong to same segment.
//            predictions.insert(category);

            std::pair<int, bool> pair(category, true); // true means positive strand
            toReturn[rowNum].insert(pair);
//        }
    }
    return toReturn;
}

int TensorflowInference::getOutputShape() {
    return output_shape;
}
int TensorflowInference::getInputShape() {
    return input_shape;
}

int TensorflowInference::getK() {
    return K;
}

void TensorflowInference::close() {
    session->Close();
}