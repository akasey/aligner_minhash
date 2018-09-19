//
// Created by Akash Shrestha on 9/14/18.
//

#ifndef ALIGNER_MINHASH_LONG_PRIORITYQUEUE_WRAPPER_H
#define ALIGNER_MINHASH_LONG_PRIORITYQUEUE_WRAPPER_H

#include "common.h"
#include "minhash.h"
#include <queue>

class LongPriorityQueueWrapper {
typedef uint32_t KeyType;
private:
    struct Element {
        float score;
        KeyType key;

        bool operator<(const Element &rhs) const {
            if (score == rhs.score)
                return key < rhs.key;
            return score < rhs.score;
        }

        bool operator==(const Element &rhs) const {
            return score == rhs.score && key==rhs.key;
        }
    };

    std::map<KeyType, std::set<Minhash::Neighbour> *> queueMapping;
    std::priority_queue<Element> queue;

    KeyType msbSetNumber; // const
    KeyType makeKey(int fragment, int partition, bool forwardStrand);
    void decodeKey(KeyType key, int *fragment, int *partition, bool *forwardStrand);

public:
    LongPriorityQueueWrapper();
    ~LongPriorityQueueWrapper();

    void addQueue(int fragment, int partition, bool forwardStrand, std::set<Minhash::Neighbour> *queue);
    bool hasNext();
    void pop(int *fragment, int *partition, bool *forwardStrand, Minhash::Neighbour * retNeighbour);
};


#endif //ALIGNER_MINHASH_LONG_PRIORITYQUEUE_WRAPPER_H
