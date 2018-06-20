//
// Created by Akash Shrestha on 6/19/18.
//

#ifndef ALIGNER_MINHASH_PRIORITYQUEUE_WRAPPER_H
#define ALIGNER_MINHASH_PRIORITYQUEUE_WRAPPER_H

#include "common.h"
#include "minhash.h"
#include <queue>

/**
 * Let's say neural network predicts the read as belonging to two segments. When searching minhash, we search both positive and negative strand of read on that segment. Total set<Minhash::Neighbour>s = 4;
 * Those four sets are sorted and first element will contain Neighbour with highest score. This class wraps a priority_queue which will select the highest score out of those four sets.
 * Priority queue stores the beginning element of each sets as struct Element. Element->key is unsigned 16bit int. MSB is set if that element comes from set of positive strand and other bits are for NN predicted segment.
 * Here the assumption is number of predicted segments will be less than 2^16. makeKey() and decodeKey() are for those bit manipulations.
 *
 * To add new set, use addQueue(). It will add first element of set to priority queue.
 * To pop an element, use pop(). It will pop max element of priority queue and remove that element from corresponding set. Also it takes the new beginning element of that set and put it into priority queue.
 */
class PriorityQueueWrapper {
typedef uint16_t KeyType;


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
    KeyType makeKey(int partition, bool forwardStrand);
    void decodeKey(KeyType key, int *partition, bool *forwardStrand);

public:
    PriorityQueueWrapper();
    ~PriorityQueueWrapper();

    void addQueue(int partition, bool forwardStrand, std::set<Minhash::Neighbour> *queue);
    bool hasNext();
    void pop(int *partition, bool *forwardStrand, Minhash::Neighbour * retNeighbour);

};


#endif //ALIGNER_MINHASH_PRIORITYQUEUE_WRAPPER_H
