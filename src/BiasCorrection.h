#ifndef BIASCORRECTION_H
#define BIASCORRECTION_H

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/algorithm/string.hpp>
#include "KmerHash.h"
#include "EMAlgorithm.h"

using namespace std;

class BiasCorrection {
    public:
        int NG;
        vector <int> NE;
        vector <int> NX;
        vector <vector<pair<int, int> > > geneBoundary;
        string genomePath;
        int readLength;

        int NC;
        vector <Eigen::MatrixXd> X;
        vector <Eigen::MatrixXd> Y;
        vector <Eigen::MatrixXd> Beta;
        vector <Eigen::MatrixXd> P;

        BiasCorrection();
        BiasCorrection(KmerHash&, EMAlgorithm&);
        ~BiasCorrection();

        inline int parseBP(char c);
        void initialX();
        void initialY(EMAlgorithm&);
        void PCA(int);
        void regression(int);
        void biasCorrection();

        void work();
};

#endif
