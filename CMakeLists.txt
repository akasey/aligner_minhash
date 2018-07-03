cmake_minimum_required(VERSION 3.8)
project(aligner_minhash)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG_MODE=1")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DDEBUG_MODE=0")

# tensorflow package
set(TensorflowCC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tf-standalone/")
find_package(TensorflowCC REQUIRED)

# SSW(Complete-Striped-Smith-Waterman-Library) Library
set(SswCC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/")
find_package(SswCC REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(include)

add_executable(minhash_index minhash_index.cc include/minhash.cc include/common.cc include/kmer.cc include/indexer_jobparser.cc include/fastamm.cc include/ThreadPool.h)
target_link_libraries(minhash_index ${CMAKE_THREAD_LIBS_INIT})

add_executable(aligner align.cc include/minhash.cc include/common.cc
        include/kmer.cc include/tf_metaparser.cc include/fastamm.cc
        include/ThreadPool.h include/tensorflow_inference.cc include/fastaQ.cc
        include/input_reader.cc include/sam_writer.cc include/indexer_jobparser.cc
        include/priorityqueue_wrapper.cc
        dependencies/ssw/src/ssw_cpp.cpp dependencies/ssw/src/ssw.c)
target_link_libraries(aligner TensorflowCC::Static)
target_link_libraries(aligner ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(aligner SswCC::Static)


add_executable(minhash_test tests/MinhashTest.cc include/minhash.cc include/common.cc include/kmer.cc
        dependencies/ssw/src/ssw_cpp.cpp dependencies/ssw/src/ssw.c)
target_link_libraries(minhash_test SswCC::Static)
target_link_libraries(minhash_test ${CMAKE_THREAD_LIBS_INIT})

add_executable(fastq_test tests/FastaQTest.cc include/common.cc include/fastaQ.cc include/input_reader.cc)

add_executable(inference_test tests/TensorflowInferenceTest.cc include/minhash.cc include/common.cc include/kmer.cc include/tensorflow_inference.cc include/tf_metaparser.cc)
target_link_libraries(inference_test TensorflowCC::Static)
target_link_libraries(inference_test ${CMAKE_THREAD_LIBS_INIT})

add_executable(playground tests/PlayGround.cc include/priorityqueue_wrapper.cc include/minhash.cc include/common.cc include/kmer.cc)
target_link_libraries(playground SswCC::Static)