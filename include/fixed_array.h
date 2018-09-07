//
// Created by Akash Shrestha on 9/6/18.
//

#ifndef ALIGNER_MINHASH_FIXED_ARRAY_H
#define ALIGNER_MINHASH_FIXED_ARRAY_H

#include "common.h"

template<typename T>
class FixedArray {
private:
    size_t definedSize;
    size_t writeHead;
    std::shared_ptr<T> elements;
public:
    FixedArray(size_t _size) {
        definedSize = _size;
        writeHead = 0;
        elements = std::shared_ptr<T> (new T[definedSize]);
    }

    ~FixedArray() {

    }

    void push_back(T elem) {
        if (writeHead < definedSize) {
            elements.get()[writeHead] = elem;
            writeHead++;
        }
        else {
            throw std::out_of_range("Size: " + std::to_string(definedSize) + " WriteHead: " + std::to_string(writeHead) );
        }
    }

    T at(size_t idx) {
        if (idx < definedSize) {
            return elements.get()[idx];
        }
        else {
            throw std::out_of_range("Size: " + std::to_string(definedSize) + " Access: " + std::to_string(idx) );
        }
    }

    size_t capacity() {
        return definedSize;
    }

    size_t size() {
        return writeHead;
    }
};


#endif //ALIGNER_MINHASH_FIXED_ARRAY_H
