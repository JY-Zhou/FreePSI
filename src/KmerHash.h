#ifndef KMERHASH_H
#define KMERHASH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <ctime>
#include <tuple>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

using namespace std;

class KmerHash {
    public:
        string genomePath;
        string exonBoundaryPath;
        vector <string> allReadPaths;

        bool strandSpec;
        double BPQualBound, AVGQualBound;
        int lowExpKmerBound;
        long long mask, shiftLeftMost;
        vector <int> chromST;
        vector <string> chromName;

        int K;
        int readLength;
        long long NW;
        int NG;
        vector <int> NE;

        unordered_map <long long, long long> kmerCount;
        unordered_map <long long, unordered_map<long long, double> > kmerTable;
        vector <vector<pair<int, int> > > geneBoundary;
        vector <string> geneStrand;

        KmerHash();
        KmerHash(int, string, string, vector <string>, int, string);
        ~KmerHash();
        
        void readReads(vector <string>);
        void readReadsFromFastq(vector <string>);
        void readReadsFromJellyfish(vector <string>);
        void readGenome(string, string);
        void readFromBed(string, string);
        void buildKmerTable(string);
        void mergeKmerTable();

        inline long long encryptExonId(int, int, int);

        inline long long parseBP(char);
        inline long long parseAntiBP(char);
        inline double positiveContribution(int, int, int, int, int);
        inline double negativeContribution(int, int, int, int, int);
        inline void insertKmerTable(long long, long long, double);
        inline void updateKmerCount(long long, long long, long long&);
};

#endif
