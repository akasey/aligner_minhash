//
// Created by Akash Shrestha on 6/11/18.
//

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include "priorityqueue_wrapper.h"

std::set<Minhash::Neighbour> makeRandom() {
    std::set<Minhash::Neighbour> neigh;
    int max = rand() % 10 + 1;
    for (int i=0; i< max; i++) {
        Minhash::Neighbour n;
        n.id = i+1;
        n.jaccard = rand() % 10;
        neigh.insert(n);
        std::cout << "NewNeighbour: (" << n.id << "," << n.jaccard << ")\n";
    }
    return neigh;
}

//	Align a pair of genome sequences.
int main (int argc, char * const argv[]) {
    PriorityQueueWrapper wrapper;
    srand(123);
    int max = 5;
    for (int i=0; i<max; i++) {
        std::set<Minhash::Neighbour> queue1 = makeRandom();
        wrapper.addQueue(i, i%2, &queue1);
    }


    while(wrapper.hasNext()) {
        Minhash::Neighbour neighbour;
        int partition;
        bool forwardStrand;
        wrapper.pop(&partition, &forwardStrand, &neighbour);
        std::cout << "partition: " << partition << " forwardStrand: " << forwardStrand << " neighbour:(" << neighbour.id << " ," << neighbour.jaccard  << ")" << std::endl;
    }
}