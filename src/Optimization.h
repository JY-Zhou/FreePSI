#ifndef OPTIMIZATION_H
#define OPTIMIZATION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdio>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "EMAlgorithm.h"

using namespace std;

extern "C" void slsqp_(int* M, int* MEQ, int* LA, int* N,
                        double* X, double* XL, double* XU, 
                        double* F, double* C, double* G, double* A,
                        double* ACC, int* ITER, int* MODE,
                        double* W, int* L_W, int* JW, int* L_JW);

class EMAlgorithm;

class Optimization {
    public:
        EMAlgorithm& emPars;
        omp_lock_t tlock;

        Optimization(EMAlgorithm&);
        ~Optimization();
        
        Eigen::MatrixXd optimizeQ(int, double, int);
        double QFunction(Eigen::MatrixXd&, int);
        Eigen::MatrixXd QGradient(Eigen::MatrixXd&, int);
        Eigen::MatrixXd QConstraints(Eigen::MatrixXd&, int);
        Eigen::MatrixXd QConstraintsNormal(Eigen::MatrixXd&, int);
};

#endif
