//
// Created by Akash Shrestha on 5/7/18.
//

#ifndef ALIGNER_MINHASH_COMMON_H
#define ALIGNER_MINHASH_COMMON_H

#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <set>
#include <memory>
#include <vector>
#include <map>
#include <exception>
#include <cstdlib>
#include <dirent.h>
#include "logger.h"

#define NULL_CHECK(val, message)  if (val == NULL) { LOG(ERROR) << message; exit(-69); }

#define PROGRAM_NAME "asm_aligner"
#define PAIRED_DISTANCE_THRESHOLD 500

struct ElementNotFoundException : public std::exception {
    std::string message;
    ElementNotFoundException(std::string msg) : message(msg) {}
    const char * what () const throw () {
        std::string m = "Element not found " + message;
        return m.c_str();
    }
};

struct TensorflowInferenceException : public std::exception {
    std::string message;
    TensorflowInferenceException(std::string msg) : message(msg) {}
    const char * what () const throw () {
        std::string m = "TFInferenceExcpt::" + message;
        return m.c_str();
    }
};

void write_in_file(const void *ptr, size_t size, size_t count, FILE *stream);

void read_from_file(void *ptr, size_t size, size_t count, FILE *stream);

std::string & trim(std::string & str);
std::string & ltrim(std::string & str);
std::string & rtrim(std::string & str);
bool str_endswith(std::string &original, std::string &with);
std::vector<std::string> split(std::string sentence, std::string delimiter);

/* File system related */
std::vector<std::string> getFilesInDirectory(std::string directory, std::string extensionFilter = "");
void mkdir(std::string dirname);
std::string basename(std::string &fqdn);

#endif //ALIGNER_MINHASH_COMMON_H
