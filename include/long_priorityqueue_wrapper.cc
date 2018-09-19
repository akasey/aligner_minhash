//
// Created by Akash Shrestha on 9/14/18.
//

#include "long_priorityqueue_wrapper.h"


LongPriorityQueueWrapper::LongPriorityQueueWrapper() : msbSetNumber(0) {
    // Make MSB = 1 for PositiveStrand because 2^16=65536.. Assuming there won't be 65534 partitions in genome.
    // Otherwise uint32_t for KeyType
    KeyType n = 0;
    n = ~n;
    KeyType x = n >> 1;
    n = ~(n & x);
    // n has MSB set to 1
    msbSetNumber = n;
}

LongPriorityQueueWrapper::~LongPriorityQueueWrapper() {
    for (std::map<KeyType, std::set<Minhash::Neighbour> *>::iterator itr=queueMapping.begin(); itr!=queueMapping.end(); itr++) {
        delete itr->second;
    }
    queueMapping.clear();
    // queue.clear;
}

LongPriorityQueueWrapper::KeyType LongPriorityQueueWrapper::makeKey(int fragment, int partition, bool forwardStrand) {
    // 1st bit forwardStrand, 15-bit for fragment,  partition 16-bit
    KeyType key = (KeyType) fragment;
    key = key << 16;
    key = key | (KeyType) partition;
    KeyType n = 0;
    if (forwardStrand) {
        n = msbSetNumber;
    }
    return key | n;
}

void LongPriorityQueueWrapper::decodeKey(KeyType key, int *fragment, int *partition, bool *forwardStrand) {
    *forwardStrand = (key & msbSetNumber);
    KeyType negmask = (~(uint16_t)0 << 16);
    *fragment = (key & (~msbSetNumber) & negmask) >> 16;
    negmask = negmask >> 16;
    *partition = key & (~msbSetNumber) & negmask;
}

void LongPriorityQueueWrapper::addQueue(int fragment, int partition, bool forwardStrand, std::set<Minhash::Neighbour> *neighourSet) {
    if (neighourSet->size() > 0) {
        KeyType key = makeKey(fragment, partition, forwardStrand);
        queueMapping[key] = new std::set<Minhash::Neighbour>();
        queueMapping[key]->insert(neighourSet->begin(), neighourSet->end());
        Element ele;
        ele.key = key;
        ele.score = neighourSet->begin()->jaccard;
//    std::cout << "PushElement: (" << ele.key << "," << ele.score << ")\n";
        queue.push(ele);
    }
}

bool LongPriorityQueueWrapper::hasNext() {
    return !queue.empty();
}

void LongPriorityQueueWrapper::pop(int *fragment, int *partition, bool *forwardStrand, Minhash::Neighbour * retNeighbour) {
    Element ele = queue.top();
//    std::cout << "PoppedElement: (" << ele.key << "," << ele.score << ")\n";
    queue.pop();
    std::set<Minhash::Neighbour> *neighourSet = queueMapping[ele.key];
    *retNeighbour = *(neighourSet->begin());
    decodeKey(ele.key, fragment, partition, forwardStrand);
    neighourSet->erase(neighourSet->begin());

    if (neighourSet->size()>0) {
        ele.score = neighourSet->begin()->jaccard;
        queue.push(ele);
    }
}