#ifndef ROSENGRADIENTDESCEND_H
#define ROSENGRADIENTDESCEND_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include <cstdio>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/PardisoSupport>

#include "EMAlgorithm.h"

using namespace std;

class EMAlgorithm;

class RosenGradientDescend{
    private:
        const int NE, NA, NX, Njunc, g;
        const EMAlgorithm& emPars;
        const double GoldenSplit = (sqrt(5.0) - 1.0) / 2.0;
        double epsZero, epsCtrl;

    public:
        static vector<double> timer;

        RosenGradientDescend(EMAlgorithm&, int);
        ~RosenGradientDescend();
        
        Eigen::MatrixXd optimizeQ(double, int);
        double lineSearch(Eigen::MatrixXd&, Eigen::MatrixXd&, double, double);
        double QFunction(Eigen::MatrixXd&);
        Eigen::MatrixXd QGradient(Eigen::MatrixXd&);
};

#endif
