//
// Created by Akash Shrestha on 7/19/18.
//

#ifndef ALIGNER_MINHASH_PAIRED_PRIORITYQUEUE_WRAPPER_H
#define ALIGNER_MINHASH_PAIRED_PRIORITYQUEUE_WRAPPER_H

#include "common.h"
#include "minhash.h"
#include "iterator"
#include <queue>

struct QueueEmptyException : public std::exception {
    std::string message;
    QueueEmptyException(std::string msg) : message(msg) {}
    const char * what () const throw () {
        std::string m = "QueueEmptyException::" + message;
        return m.c_str();
    }
};

class PairedPriorityQueueWrapper {
typedef uint32_t KeyType;

private:
    struct Element {
        float score;
        float distance; // distance
        KeyType key;

        bool operator <(const Element &rhs) const {
            if (score == rhs.score)
                return distance < rhs.distance;
            else if (score == rhs.score && distance == rhs.distance)
                return key < rhs.key;
            return score < rhs.score;
        }

        bool operator==(const Element &rhs) const {
            return score == rhs.score && distance==rhs.distance && key==rhs.key;
        }
    };

    struct MinhashNeighbourPairedSet {
        MinhashNeighbourPairedSet(std::set<Minhash::Neighbour> *, std::set<Minhash::Neighbour> *);
        ~MinhashNeighbourPairedSet();
        std::set<Minhash::Neighbour> *first, *second;
        std::set<Minhash::Neighbour>::iterator firstItr, secondItr; // one way traversal
        uint16_t firstInc, secondInc; // other way traversal.. integers to make sure we don't have
                                      // same inc in first and second coz the iterators will get those and will make pairs redundant
        void oneWayTraversalPeek(float *score, int *distance);
        void otherWayTraversalPeek(float *score, int *distance);
        std::pair<Minhash::Neighbour, Minhash::Neighbour> oneWayTraversalMove();
        std::pair<Minhash::Neighbour, Minhash::Neighbour> otherWayTraversalMove();
        bool isEmpty();
        std::pair<Minhash::Neighbour, Minhash::Neighbour> top();
        Element peek();
    };

    std::map<KeyType, MinhashNeighbourPairedSet *> queueMapping;
    std::priority_queue<Element> queue;

    KeyType firstMsbSetNumber; // const
    KeyType secondMsbSetNumber; // const
    KeyType makeKey(int firstPartition, bool firstForwardStrand, int secondPartition, bool secondForwardStrand);
    void decodeKey(KeyType key, int *firstPartition, bool *firstForwardStrand, int *secondPartition, bool *secondForwardStrand);

public:
    PairedPriorityQueueWrapper();
    ~PairedPriorityQueueWrapper();

    void addQueue(int firstPartition, bool firstForwardStrand, std::set<Minhash::Neighbour> *firstQueue,
                  int secondPartition, bool secondForwardStrand, std::set<Minhash::Neighbour> *secondQueue );
    bool hasNext();
    void pop(int *firstPartition, bool *firstForwardStrand, Minhash::Neighbour *firstNeighbour,
             int *secondPartition, bool *secondForwardStrand, Minhash::Neighbour *secondNeighbour);
};


#endif //ALIGNER_MINHASH_PAIRED_PRIORITYQUEUE_WRAPPER_H
