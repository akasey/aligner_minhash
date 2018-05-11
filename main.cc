#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include "include/minhash.h"

using namespace std;

struct Reads {
    string sequence;
    int position;
    vector<int> predicted_segments;
    vector<int> predicted_positions;
};

vector<string> split(string sentence, string delimiter){
    vector<string> list;
    size_t pos = 0;
    string token;
    while ((pos = sentence.find(delimiter)) != string::npos) {
        token = sentence.substr(0, pos);
        list.push_back(token);
        sentence.erase(0, pos + delimiter.length());
    }
    list.push_back(sentence);
    return list;
}

string trim(const string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

vector<Reads> readPredictions(string filename) {
    ifstream infile(filename);
    string line;
    vector<string> fasta;
    std::string::size_type sz;
    vector<Reads> reads;
    while(getline(infile, line)) {
        vector<string> tokens = split(line, "$$");
        Reads read;
        for (int i=0; i<tokens.size(); i++) {
            switch (i) {
                case 0: read.sequence = tokens[i]; break;
                case 1: read.position = stoi(tokens[i], &sz); break;
                case 2:
                    if (trim(tokens[i]).length() != 0) {
                        vector<string> splits = split(tokens[i], ",");
                        for (int j = 0; j < splits.size(); j++) {
                            int segment = stoi(splits[j], &sz);
                            read.predicted_segments.push_back(segment);
                        }
                    }
                    break;
            }
        }
        reads.push_back(read);
    }
    return reads;
}

void readGenome(string filename, string &genome) {
    ifstream infile(filename);
    string line;
    genome = "";
    getline(infile, line);
    while(getline(infile, line)) {
        genome += line;
    }
}

map<int, string> getSegments(string genome, int segmentLength, int windowLen, int **&coordMapping) {
    int genomeLength = genome.length();
    double numSegments = ceil((double)genomeLength / (double)segmentLength);
    coordMapping = new int*[(int)numSegments];
    map<int, string> segments;
    for (int i=0; i<numSegments; i++) {
        int start = fmax(0, i * segmentLength - floor(windowLen / 2.0));
        int end = fmin((i + 1) * segmentLength + ceil(windowLen / 2.0) - 1, genomeLength);
        coordMapping[i] = new int[2];
        coordMapping[i][0] = start; coordMapping[i][1] = end;
        string segment = genome.substr(start, end-start+1);
        segments[i] = segment;
    }
    return segments;
}

map<int, string> makeSlidingWindow(string segment, int *coord, int windowLength, int strides) {
    map<int, string> windows;
    int end_coord = segment.length() - windowLength;
    for (int start=0; ; start = min(start+strides, end_coord) ) {
        string sub = segment.substr(start, windowLength);
        windows[coord[0]+start] = sub;
        if (start == end_coord)
            break;
    }
    return windows;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <refGenome.fasta> <...predictions...>" << endl;
        exit(1);
    }

    string refGenome = argv[1];
    int segmentLength = 5000;
    int windowLength = 200;
    int K = 7;
    int strides = K-1;
    typedef MinHash<int, string> HashStore;

    string genome = "";
    cout << "Reading genome..." << endl;
    readGenome(refGenome, genome);
    int **coordMapping;
    cout << "Making segments..." << endl;
    map<int, string> segments = getSegments(genome, segmentLength, windowLength, coordMapping);
    int numSegment = segments.size();

    cout << "Minhashing.. numSegments: " << numSegment << "..." << endl;
    HashStore minhashes[numSegment];
    for (map<int,string>::iterator itr = segments.begin(); itr != segments.end(); itr++) {
        map<int, string> windows = makeSlidingWindow(itr->second, coordMapping[itr->first], windowLength, strides);
        for (map<int, string>::iterator itrr = windows.begin(); itrr != windows.end(); itrr++) {
            minhashes[itr->first].addDocument(itrr->first, itrr->second);
        }
    }

    for (int argumentIdx=2; argumentIdx < argc; argumentIdx++) {
        string predictionFile = argv[argumentIdx];
        std::string baseDir = predictionFile.substr(0, predictionFile.find_last_of('/'));
        cout << "Predictions File: " << predictionFile << endl;
        cout << "Output File: " << baseDir + "/predictedLocations" << endl;
        cout << "Predicting..." << endl;
        ofstream outfile(baseDir + "/predictedLocations");
        vector<Reads> reads = readPredictions(predictionFile);
        for (int i = 0; i < reads.size(); i++) {
            Reads read = reads[i];
            set<HashStore::Neighbour> allNeighbours;
            for (vector<int>::iterator itr = read.predicted_segments.begin();
                 itr != read.predicted_segments.end(); itr++) {
                int segment = *itr;
                if (segment < numSegment) {
                    set<HashStore::Neighbour> neighbours = minhashes[segment].findNeighbours(read.sequence);
                    allNeighbours.insert(neighbours.begin(), neighbours.end());
                }
            }
            outfile << read.sequence << "$$" << read.position << "$$";
            for (set<HashStore::Neighbour>::iterator itr = allNeighbours.begin(); itr != allNeighbours.end(); itr++) {
                outfile << itr->id << "(" << itr->jaccard << "),";
            }
            outfile << endl;
        }
        outfile.close();
    }

    return 0;
}