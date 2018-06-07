cmake_minimum_required(VERSION 3.8)
project(aligner_minhash)

set(CMAKE_CXX_STANDARD 11)
set(TensorflowCC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tf-standalone/")
find_package(TensorflowCC REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(include)

add_executable(minhash_index minhash_index.cc include/minhash.cc include/common.cc include/kmer.cc include/indexer_jobparser.cc include/fastamm.cc include/ThreadPool.h)
target_link_libraries(minhash_index ${CMAKE_THREAD_LIBS_INIT})

add_executable(aligner align.cc include/minhash.cc include/common.cc include/kmer.cc include/tf_metaparser.cc include/fastamm.cc include/ThreadPool.h include/tensorflow_inference.cc include/fastaQ.cc include/input_reader.cc)
target_link_libraries(aligner TensorflowCC::Static)
target_link_libraries(aligner ${CMAKE_THREAD_LIBS_INIT})


add_executable(minhash_test tests/MinhashTest.cc include/minhash.cc include/common.cc include/kmer.cc)
add_executable(fastq_test tests/FastaQTest.cc include/common.cc include/fastaQ.cc include/input_reader.cc)

add_executable(inference_test tests/TensorflowInferenceTest.cc include/minhash.cc include/common.cc include/kmer.cc include/tensorflow_inference.cc include/tf_metaparser.cc)
target_link_libraries(inference_test TensorflowCC::Static)
target_link_libraries(inference_test ${CMAKE_THREAD_LIBS_INIT})