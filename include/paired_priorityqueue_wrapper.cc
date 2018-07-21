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
    return firstItr == first->end() && secondInc == second->size();
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

std::pair<Minhash::Neighbour, Minhash::Neighbour> PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::oneWayTraversalMove() {
    std::pair<Minhash::Neighbour, Minhash::Neighbour> toRet(*firstItr, *secondItr);
    do {
        secondItr++;
        if (secondItr == second->end()) {
            firstItr++;
            secondItr = second->begin();
        }
        if (firstItr == first->end())
            throw new QueueEmptyException("OneWayMove: I'm out");
    } while ( ! abs(firstItr->id - secondItr->id) <= PAIRED_DISTANCE_THRESHOLD );
    return toRet;
}

void PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::otherWayTraversalMove() {
    std::pair<Minhash::Neighbour, Minhash::Neighbour> toRet(*(first->begin()+firstInc), *(second->begin()+secondInc));
    do {
        firstInc++;
        if (firstInc == first->size()) {
            secondInc++;
            firstInc = 0;
        }
        if (secondInc == second->size())
            throw new() QueueEmptyException("OtherWayMove: I'm out");
    } while ( firstInc != second && ! abs((first->begin()+firstInc)->id - (second->begin()+secondInc)->id) <= PAIRED_DISTANCE_THRESHOLD );
    return toRet;
}

std::pair<Minhash::Neighbour, Minhash::Neighbour> PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::top() {
    float oneWayScore; int oneWayDistance;
    float otherWayScore; int otherWayDistance;
    oneWayTraversalPeek(&oneWayScore, &oneWayDistance);
    otherWayTraversalPeek(&otherWayScore, &otherWayDistance);
    if (oneWayScore < otherWayScore)
        return oneWayTraversalMove();
    else if (otherWayScore > oneWayScore)
        return otherWayTraversalMove();
    else if (oneWayScore == otherWayScore && oneWayDistance < otherWayDistance)
        return oneWayTraversalMove();
    else
        return otherWayTraversalMove();
}


PairedPriorityQueueWrapper::PairedPriorityQueueWrapper() {
    KeyType n= 0;
    n = ~n;
    uint16_t x = n >> 1;
    n = ~(n & x);
    firstMsbSetNumber = n;
    n = n >> 1;
    secondMsbSetNumber = n;
}

PairedPriorityQueueWrapper::~PairedPriorityQueueWrapper() {
    queueMapping.clear();
}

KeyType PairedPriorityQueueWrapper::makeKey(int firstPartition, bool firstForwardStrand, int secondPartition, bool secondForwardStrand) {
    KeyType n = firstPartition;
    n = n << 15 | secondPartition;
    if (firstForwardStrand)
        n = n | firstMsbSetNumber;
    if (secondForwardStrand)
        n = n | secondMsbSetNumber;
    return n;
}

void PairedPriorityQueueWrapper::decodeKey(KeyType key, int *firstPartition, bool *firstForwardStrand, int *secondPartition, bool *secondForwardStrand) {
    firstPartition = ( key & (~firstMsbSetNumber) & (~secondMsbSetNumber) ) >> 15;
    firstForwardStrand = key & firstMsbSetNumber;

    secondPartition = ( key & (~firstMsbSetNumber) & (~secondMsbSetNumber) )
}