//
// Created by Akash Shrestha on 9/6/18.
//

#include <iostream>
#include "fixed_array.h"

int main() {
    FixedArray<int> fixedArray(10);

    for (int i=0; i<fixedArray.capacity(); i++) {
        fixedArray.push_back(i*i);
    }

    for (int i=0; i<fixedArray.size(); i++) {
        std::cout << fixedArray.at(i) << std::endl;
    }
    return 0;
}