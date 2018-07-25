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
        throw QueueEmptyException("OneWay: I'm out");
    *score = firstItr->jaccard + secondItr->jaccard;
    *distance = abs((int)firstItr->id - (int)secondItr->id);
}

void PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::otherWayTraversalPeek(float *score, int *distance) {
    if (firstInc == first->size())
        throw QueueEmptyException("OtherWay: I'm out");
    auto firstNav = std::next(first->begin(), firstInc);
    auto secondNav = std::next(second->begin(), secondInc);
    *score = firstNav->jaccard + secondNav->jaccard;
    *distance = abs((int)firstNav->id - (int)secondNav->id);
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
            throw QueueEmptyException("OneWayMove: I'm out");
    } while ( ! abs((int)firstItr->id - (int)secondItr->id) <= PAIRED_DISTANCE_THRESHOLD );
    return toRet;
}

std::pair<Minhash::Neighbour, Minhash::Neighbour> PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::otherWayTraversalMove() {
    auto firstNav = std::next(first->begin(), firstInc);
    auto secondNav = std::next(second->begin(), secondInc);
    std::pair<Minhash::Neighbour, Minhash::Neighbour> toRet(*firstNav, *secondNav);
    do {
        firstInc++;
        if (firstInc == first->size()) {
            secondInc++;
            firstInc = 0;
        }
        if (secondInc == second->size())
            throw QueueEmptyException("OtherWayMove: I'm out");
        firstNav = std::next(first->begin(), firstInc);
        secondNav = std::next(second->begin(), secondInc);
    } while (! abs((int)firstNav->id - (int)secondNav->id) <= PAIRED_DISTANCE_THRESHOLD );
    return toRet;
}

std::pair<Minhash::Neighbour, Minhash::Neighbour> PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::top() {
    float oneWayScore; int oneWayDistance;
    float otherWayScore; int otherWayDistance;
    oneWayTraversalPeek(&oneWayScore, &oneWayDistance);
    otherWayTraversalPeek(&otherWayScore, &otherWayDistance);
    Element element1, element2;
    element1.key = 1; element1.score = oneWayScore; element1.distance = oneWayDistance;
    element2.key = 1; element2.score = oneWayScore; element2.distance = oneWayDistance;
    if (element2 < element1)
        return oneWayTraversalMove();
    else
        return otherWayTraversalMove();
}

PairedPriorityQueueWrapper::Element PairedPriorityQueueWrapper::MinhashNeighbourPairedSet::peek() {
    float oneWayScore; int oneWayDistance;
    float otherWayScore; int otherWayDistance;
    oneWayTraversalPeek(&oneWayScore, &oneWayDistance);
    otherWayTraversalPeek(&otherWayScore, &otherWayDistance);
    Element element1, element2;
    element1.key = 0; element1.score = oneWayScore; element1.distance = oneWayDistance;
    element2.key = 0; element2.score = oneWayScore; element2.distance = oneWayDistance;
    if (element2 < element1)
        return element1;
    else
        return element2;
}


PairedPriorityQueueWrapper::PairedPriorityQueueWrapper() {
    KeyType n= 0;
    n = ~n;
    KeyType x = n >> 1;
    n = ~(n & x);
    firstMsbSetNumber = n;
    n = n >> 1;
    secondMsbSetNumber = n;
}

PairedPriorityQueueWrapper::~PairedPriorityQueueWrapper() {
    queueMapping.clear();
}

PairedPriorityQueueWrapper::KeyType PairedPriorityQueueWrapper::makeKey(int firstPartition, bool firstForwardStrand, int secondPartition, bool secondForwardStrand) {
    KeyType n = (uint16_t) firstPartition;
    n = n << 15 | (uint16_t) secondPartition;
    if (firstForwardStrand)
        n = n | firstMsbSetNumber;
    if (secondForwardStrand)
        n = n | secondMsbSetNumber;
    return n;
}

void PairedPriorityQueueWrapper::decodeKey(KeyType key, int *firstPartition, bool *firstForwardStrand, int *secondPartition, bool *secondForwardStrand) {
    KeyType negmask = (~(uint16_t)0 << 15);
    KeyType mask = negmask & (~firstMsbSetNumber) & (~secondMsbSetNumber);
    *firstPartition = (key & mask) >> 15;
    *firstForwardStrand = key & firstMsbSetNumber;

    negmask = ~(~0 << 15);
    *secondPartition = (key & negmask);
    *secondForwardStrand = key & secondMsbSetNumber;
}

void PairedPriorityQueueWrapper::addQueue(int firstPartition, bool firstForwardStrand, std::set<Minhash::Neighbour> *firstQueue,
              int secondPartition, bool secondForwardStrand, std::set<Minhash::Neighbour> *secondQueue ) {
    if (firstQueue->size() > 0 && secondQueue->size() > 0) {
        KeyType key = makeKey(firstPartition, firstForwardStrand, secondPartition, secondForwardStrand);
        PairedPriorityQueueWrapper::MinhashNeighbourPairedSet *pairedSet = new MinhashNeighbourPairedSet(firstQueue, secondQueue);
        queueMapping[key] = pairedSet;

        try {
            Element ele = pairedSet->peek();
            ele.key = key;
            queue.push(ele);
        }
        catch (QueueEmptyException ex) {

        }
    }
}

bool PairedPriorityQueueWrapper::hasNext() {
    return !queue.empty();
}

void PairedPriorityQueueWrapper::pop(int *firstPartition, bool *firstForwardStrand, Minhash::Neighbour *firstNeighbour,
         int *secondPartition, bool *secondForwardStrand, Minhash::Neighbour *secondNeighbour) {
    Element ele = queue.top();
    queue.pop();

    try {
        MinhashNeighbourPairedSet *pairedSet = queueMapping[ele.key];
        decodeKey(ele.key, firstPartition, firstForwardStrand, secondPartition, secondForwardStrand);
        std::pair<Minhash::Neighbour, Minhash::Neighbour> pairedNeighbours = pairedSet->top();
        *firstNeighbour = pairedNeighbours.first;
        *secondNeighbour = pairedNeighbours.second;

        Element element = pairedSet->peek();
        element.key = ele.key;
        queue.push(element);
    }
    catch (QueueEmptyException ex) {

    }

}