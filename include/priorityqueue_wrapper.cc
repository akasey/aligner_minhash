//
// Created by Akash Shrestha on 6/19/18.
//

#include "priorityqueue_wrapper.h"

PriorityQueueWrapper::PriorityQueueWrapper() : msbSetNumber(0) {
    // Make MSB = 1 for PositiveStrand because 2^16=65536.. Assuming there won't be 65534 partitions in genome.
    // Otherwise uint32_t for KeyType
    KeyType n = 0;
    n = ~n;
    KeyType x = n >> 1;
    n = ~(n & x);
    // n has MSB set to 1
    msbSetNumber = n;
}

PriorityQueueWrapper::~PriorityQueueWrapper() {
    for (std::map<KeyType, std::set<Minhash::Neighbour> *>::iterator itr=queueMapping.begin(); itr!=queueMapping.end(); itr++) {
        delete itr->second;
    }
    queueMapping.clear();
    // queue.clear;
}


PriorityQueueWrapper::KeyType PriorityQueueWrapper::makeKey(int partition, bool forwardStrand) {
    KeyType key = (KeyType) partition;
    KeyType n = 0;
    if (forwardStrand) {
        n = msbSetNumber;
    }
    return key | n;
}

void PriorityQueueWrapper::decodeKey(KeyType key, int *partition, bool *forwardStrand) {
    *forwardStrand = (key & msbSetNumber);
    *partition = key & (~msbSetNumber);
}

void PriorityQueueWrapper::addQueue(int partition, bool forwardStrand, std::set<Minhash::Neighbour> *neighourSet) {
    if (neighourSet->size() > 0) {
        KeyType key = makeKey(partition, forwardStrand);
        queueMapping[key] = new std::set<Minhash::Neighbour>();
        queueMapping[key]->insert(neighourSet->begin(), neighourSet->end());
        Element ele;
        ele.key = key;
        ele.score = neighourSet->begin()->jaccard;
//    std::cout << "PushElement: (" << ele.key << "," << ele.score << ")\n";
        queue.push(ele);
    }
}

bool PriorityQueueWrapper::hasNext() {
    return !queue.empty();
}

void PriorityQueueWrapper::pop(int *partition, bool *forwardStrand, Minhash::Neighbour * retNeighbour) {
    Element ele = queue.top();
//    std::cout << "PoppedElement: (" << ele.key << "," << ele.score << ")\n";
    queue.pop();
    std::set<Minhash::Neighbour> *neighourSet = queueMapping[ele.key];
    *retNeighbour = *(neighourSet->begin());
    decodeKey(ele.key, partition, forwardStrand);
    neighourSet->erase(neighourSet->begin());

    if (neighourSet->size()>0) {
        ele.score = neighourSet->begin()->jaccard;
        queue.push(ele);
    }
}