//
// Created by Akash Shrestha on 7/19/18.
//

#ifndef ALIGNER_MINHASH_PAIRED_PRIORITYQUEUE_WRAPPER_H
#define ALIGNER_MINHASH_PAIRED_PRIORITYQUEUE_WRAPPER_H

#include "common.h"
#include "minhash.h"

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
    struct MinhashNeighbourPairedSet {
        MinhashNeighbourPairedSet(std::set<Minhash::Neighbour> *, std::set<Minhash::Neighbour> *);
        ~MinhashNeighbourPairedSet();
        std::set<Minhash::Neighbour> *first, *second;
        std::set<Minhash::Neighbour>::iterator firstItr, secondItr; // one way traversal
        uint16_t firstInc, secondInc; // other way traversal.. integers to make sure we don't have
                                      // same inc in first and second coz the iterators will get those and will make pairs redundant
        void oneWayTraversalPeek(float *score, int *distance);
        void otherWayTraversalPeek(float *score, int *distance);
        void oneWayTraversalMove(float *score, int *distance);
        void otherWayTraversalMove(float *score, int *distance);
        bool isEmpty();
        std::pair<Minhash::Neighbour, Minhash::Neighbour> top();
    };

    struct Element {
        float score;
        float distance; // distance
        KeyType key;

        bool operator <(const Element &rhs) const {
            if (score == rhs.score)
        }
    };
};


#endif //ALIGNER_MINHASH_PAIRED_PRIORITYQUEUE_WRAPPER_H
