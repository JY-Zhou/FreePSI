#ifndef CONJUGATEGRADIENTPROJECTION_H
#define CONJUGATEGRADIENTPROJECTION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <stack>

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

#include "EMAlgorithm.h"

using namespace std;

class EMAlgorithm;

class ConjugateGradientProjection{
    private:
        const int NX, NEQ, g;
        const double GoldenSplit = (sqrt(5.0) - 1.0) / 2.0;
        const double epsZero;

        EMAlgorithm & emPars;

        Eigen::MatrixXd eq;
        Eigen::SparseMatrix<double> ieq;

        Eigen::MatrixXd INV, P, H, X, grad;
        Eigen::SparseMatrix<double> N, CC, CM;

        vector <int> stateOfIeqs; //Start from 0, -1 indicate not active
        //vector <bool> dropped; //Start from 0, -1 indicate not active
        vector <int> hyperplaneSeq; //First NEQ are equalities, with value -1, other are idx of inequalities,

    public:
        static vector<double> timer, exitPort;
        static double QFunc_t, QGrad_t, ConsAdd_t;
        static int QFunc_n, QGrad_n, ConsAdd_n, iteration;

        ConjugateGradientProjection(EMAlgorithm&, int);
        ~ConjugateGradientProjection();

        Eigen::MatrixXd optimizeQ(double, int);

        void initOptimization();
        void swapLastHyperplane(int);
        void dropLastHyperplane();
        bool addHyperplane(int);
        //void updateH(Eigen::SparseMatrix<double> &);
        Eigen::SparseMatrix<double> correctInfeasible(Eigen::SparseMatrix<double> &);

        double GoldenSelection(Eigen::SparseMatrix<double>&, Eigen::SparseMatrix<double>&, double, double);
        double SecantMethod(Eigen::MatrixXd &, Eigen::MatrixXd &, double, double);
        double QFunction(Eigen::SparseMatrix<double>&);
        Eigen::MatrixXd QGradient(Eigen::MatrixXd &);
};

#endif
