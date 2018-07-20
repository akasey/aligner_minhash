//
// Created by Akash Shrestha on 7/19/18.
//

#include "paired_priorityqueue_wrapper.h"

PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::MinhashNeighbourPairedSet(std::set<Minhash::Neighbour> * _first,
                                                                                 std::set<Minhash::Neighbour> * _second) {
    first = new std::set<Minhash::Neighbour>();
    first->insert(_first->begin(), _first->end());
    second = new std::set<Minhash::Neighbour>();
    second->insert(_second->begin(), _second->end());

    firstItr = first->begin(), secondItr = second->begin(); // one way traversal
    secondInc = 0, firstInc = 1; // other way traversal
}

PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::~MinhashNeighbourPairedSet() {
    delete first;
    delete second;
}

bool PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::isEmpty() {
    return firstItr == first->end();
}

void PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::oneWayTraversalPeek(float *score, int *distance) {
    if (firstItr == first->end())
        throw new QueueEmptyException("OneWay: I'm out");
    *score = firstItr->jaccard + secondItr->jaccard;
    *distance = abs(firstItr->id - secondItr->id);
}

void PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::otherWayTraversalPeek(float *score, int *distance) {
    if (firstInc == first->size())
        throw new QueueEmptyException("OtherWay: I'm out");
    *score = (first->begin()+firstInc)->jaccard + (second->begin()+secondInc)->jaccard;
    *distance = abs((first->begin()+firstInc)->id - (second->begin()+secondInc)->id);
}

std::pair<Minhash::Neighbour, Minhash::Neighbour> PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::top() {
    while ( ! fabs(firstItr->id - secondItr->id) <= PAIRED_DISTANCE_THRESHOLD ) {
        if (isEmpty())
            throw new QueueEmptyException("I'm out");
        secondItr++;
        if (secondItr == second->end()) {
            firstItr++;
            secondItr = second->begin();
        }
    }
    return std::pair<Minhash::Neighbour, Minhash::Neighbour>(*firstItr, *secondItr);
}