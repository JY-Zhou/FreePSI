#ifndef EMALGORITHM_H
#define EMALGORITHM_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <tuple>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>

#include <omp.h>

#include "ConjugateGradientProjection.h"
#include "KmerHash.h"

using namespace std;

class ConjugateGradientProjection;

class EMAlgorithm {
    public:
        const double EPS = 1e-9;
        const int NSEGMENT = 50;
        int nThread = 16;
        int procSegment = 0;
        int exonNumBound = 40;

        KmerHash& kmerHasher;

        double lambda;
        int K;
        int readLength;

        long long NW;
        int NG;
        vector <int> NE;
        vector <int> NX;
        vector <int> NA;

        vector <Eigen::MatrixXd> L;
        vector <Eigen::SparseMatrix <double> > A;
        vector <Eigen::SparseMatrix <double> > C;
        Eigen::MatrixXd W;

        Eigen::SparseMatrix <double> Z;
        vector <int> isAvaliableZ; //0 is empty, 1 is ok, 2 contains over 'exonNumBound' exons
        vector <Eigen::SparseMatrix <double> > X;
        vector <Eigen::SparseMatrix <double> > Mu;
        double preLikelihood, newLikelihood;

        vector <ConjugateGradientProjection> CGPAOptmizer;

        vector <vector <double> > PSI;

        EMAlgorithm(KmerHash&, int nThread);
        ~EMAlgorithm();

        inline tuple <int, int, int> decryptExonId(long long);

        void initialIndices();
        void initialCoefficients(KmerHash&);
        void initialConstraints();
        void initialVariables();
        void correctInfeasibleX(int);
        void correctInfeasibleVariables();

        double computeLogLikelihood();

        void eStep();
        void mStep(double, int);
        void work(string);

        void computePSI(string);
};

#endif
