#include "Optimization.h"
#define EIGEN_USE_MKL_ALL

Optimization::Optimization(EMAlgorithm& _emPars): emPars(_emPars) {
    //omp_init_lock(&tlock);
}

Optimization::~Optimization() {
    //omp_destroy_lock(&tlock);
}

Eigen::MatrixXd Optimization::optimizeQ(int g, double eps, int iter) {
    //cout << "\nGene-" << g << " started ...." << omp_get_thread_num() << endl;
    //omp_set_lock(&tlock);
    int f_M = emPars.NA[g] + 1;
    int f_MEQ = 1;
    int f_LA = f_M;
    int f_N = emPars.NX[g];
    //Eigen::MatrixXd f_X = emPars.initialX(g);
    Eigen::MatrixXd f_X = emPars.X[g];
    Eigen::MatrixXd f_XL;
    Eigen::MatrixXd f_XU;
    f_XL.setZero(f_N, 1);
    f_XU.setOnes(f_N, 1);
    double f_F = QFunction(f_X, g);
    Eigen::MatrixXd f_C = QConstraints(f_X, g);
    Eigen::MatrixXd f_G = QGradient(f_X, g);
    Eigen::MatrixXd f_A = QConstraintsNormal(f_X, g);
    //double f_ACC = std::numeric_limits<double>::epsilon();
    double f_ACC = eps;
    int f_ITER = iter;
    int f_MODE = 0;
    int f_N1 = f_N + 1;
    int f_MINEQ = f_M - f_MEQ + 2 * f_N1;
    int f_L_W = (3 * f_N1 + f_M) * (f_N1 + 1)
              + (f_N1 - f_MEQ + 1) * (f_MINEQ + 2) + 2 * f_MINEQ
              + (f_N1 + f_MINEQ) * (f_N1 - f_MEQ) + 2 * f_MEQ + f_N1
              + f_N1 * f_N / 2 + 2 * f_M + 3 * f_N + 3 * f_N1 + 1;
    Eigen::MatrixXd f_W(f_L_W, 1);
    f_W.setZero(f_L_W, 1);
    int f_L_JW = f_MINEQ;
    Eigen::MatrixXi f_JW(f_L_JW, 1);
    f_JW.setZero(f_L_JW, 1);
    int t = 0;
    
    //ostringstream log;

    //cout << "Now at gene " << g << endl;
    //#pragma omp private (f_X, f_MODE, f_W, f_JW)
    //
    double outtime = 0;
    double intime = 0;
    double tottime = clock();

    while(true) {
        t ++;
        //if(t > 30 && emPars.NE[g] > 30) break;
        //if(t % 10 == 0) cout << "Stuck" << t << endl;
        if(std::isnan(f_X.sum())) {
            cout << endl << g << ": " <<  "-Monitor::::"<< t << ":" <<  f_X.sum() << endl;
            getchar();
        }

        //log << g << ": " <<  "-Monitor::::"<< t << ": " <<  f_X.sum() << endl;
        //log << g << ": " <<  "==" << endl;
        //log << g << ": " <<  "MODE " << f_MODE << endl;
        if(f_MODE == 0 || f_MODE == 1) {
            clock_t tbegin = clock();
            f_F = QFunction(f_X, g);
            double elap = (double)(clock() - tbegin) / CLOCKS_PER_SEC;
            if(elap > 0) {
                //cout << "QF elapse " << elap << endl;
                outtime += elap;
            }

            tbegin = clock();
            f_C = QConstraints(f_X, g);
            elap = (double)(clock() - tbegin) / CLOCKS_PER_SEC;
            if(elap > 0) {
                outtime += elap;
               // cout << "QC elapse " << elap << endl;
            }
        }

        if(f_MODE == 0 || f_MODE == -1) {
            clock_t tbegin = clock();
            f_G = QGradient(f_X, g);
            double elap = (double)(clock() - tbegin) / CLOCKS_PER_SEC;
            if(elap > 0) {
                outtime += elap;
                //cout << "QG elapse" << elap << endl;
            }

            tbegin = clock();
            f_A = QConstraintsNormal(f_X, g);
            elap = (double)(clock() - tbegin) / CLOCKS_PER_SEC;
            if(elap > 0) {
                outtime += elap;
                //cout << "QCN elapse" << elap << endl;
            }
        }

        //log << g << ": " <<  "X==" << endl;
        //log << g << ": " <<  f_X << endl;
        //log << g << ": " <<  "F==" << endl;
        //log << g << ": " <<  f_F << endl;
        //log << g << ": " <<  "C==" << endl;
        //log << g << ": " <<  f_C << endl;
        //log << g << ": " <<  "G==" << endl;
        //log << g << ": " <<  f_G << endl;
        //log << g << ": " <<  "A==" << endl;
        //log << g << ": " <<  f_A << endl;

       // printf("Before SLSQP, g = %d, F_X = %lf, [%d]\n", g, f_X.sum(), omp_get_thread_num());
       // printf("Before SLSQP, g = %d, F_G = %lf, [%d]\n", g, f_G.sum(), omp_get_thread_num());
        clock_t slsBegin = clock();
        slsqp_(&f_M, &f_MEQ, &f_LA, &f_N,
              f_X.data(), f_XL.data(), f_XU.data(),
              &f_F, f_C.data(), f_G.data(), f_A.data(),
              &f_ACC, &f_ITER, &f_MODE,
              f_W.data(), &f_L_W, f_JW.data(), &f_L_JW);
        double slsElap = (double)(clock() - slsBegin) / CLOCKS_PER_SEC;
        intime += slsElap;
        //cout << "SLSQP elapse " << slsElap << endl;
        

       // printf("After SLSQP, g = %d, F_X = %lf, [%d]\n", g, f_X.sum(), omp_get_thread_num());
       // printf("After SLSQP, g = %d, F_G = %lf, [%d]\n", g, f_G.sum(), omp_get_thread_num());

        if(f_MODE != -1 && f_MODE != 1) {
            //log << g << ": " <<  f_MODE << endl;
            //cout << emPars.NX[g] << endl;
            //cout << "Out-time is " << outtime << '\t';
            //cout << "SLSQP-time is " << intime << endl;
            break;
        }
    }
    //cout << "iter time " << t << endl;
    //cout << "x sumation " << f_X.sum() << endl;
    //
    //cout << "Optimal is " << f_F << endl;
    //cout << "Exit Mode is " << f_MODE << endl;
    //if(std::isnan(f_F) or std::isnan(f_X.sum())) {
    //    cerr << log.str() << endl;
    //} else {
    //    log.str("");
    //}

    //omp_unset_lock(&tlock);
    //cout << "TOT " << (double)(clock() - tottime) / CLOCKS_PER_SEC << endl;
    printf("%.6f\n", QFunction(f_X, g));
    return f_X;
}

double Optimization::QFunction(Eigen::MatrixXd& X, int g) {
    //Eigen::MatrixXd temp = C[g] * X;
    //temp = temp.array().log();
    //for(int w = 0; w < NW; w ++) {
    //    if(std::isinf(-temp(w, 0)) || std::isinf(temp(w, 0))) {
    //        temp(w, 0) = 0;
    //    }
    //    if(std::isnan(temp(w, 0))) {
    //        temp(w, 0) = 0;
    //    }
    //}
    //cout << "Now compute Q function" << endl;

    //double ret = - (Mu[g].transpose() * temp)(0, 0);
    //return ret;
    Eigen::SparseMatrix<double> sX = X.sparseView();
    if (omp_get_thread_num() == 1) {
        //cout << "g=" << g << " F_X=" << X << endl;
        //cout << "sX=" << sX << endl;
    }
    sX.makeCompressed();
    Eigen::SparseMatrix<double> temp = (emPars.C[g] * sX);
    //temp.makeCompressed();
    //temp = temp.array().log();
    int nz = temp.nonZeros();
    for(auto i = temp.valuePtr(); i != temp.valuePtr() + nz; ++ i) {
        if(std::isnan(*i)) *i = 0;
    }
    temp.makeCompressed();
    nz = temp.nonZeros();
    int tt = 0;
    //temp = temp.array().log();
    for(auto i = temp.valuePtr(); i != temp.valuePtr() + nz; ++ i) {
        if(*i <= 0) {
            //cout << temp.rows() << " " << temp.cols() << endl;
            //printf("heheda(0), thread number [%d], g = %d\n", omp_get_thread_num(), g);
            //cout << "SX = " << endl << sX << endl;
            //cout << "temp = " << endl << temp << endl;
            //cout << temp.sum() << endl;
            //cout << *i << endl;
            //cout << tt << endl;
            //tt ++;
        //cout << "F_X=" << X << endl;
        //cout << "sX=" << sX << endl;
            //getchar();
        }
        if(std::isnan(*i)) {
            //printf("heheda(nan), thread number [%d], g = %d\n", omp_get_thread_num(), g);
            //getchar();
        }
        if(std::isinf(*i)) {
            //printf("heheda(inf), thread number [%d], g = %d\n", omp_get_thread_num(), g);
            //getchar();
        }
        *i = log(*i);
    }
    for(auto i = temp.valuePtr(); i != temp.valuePtr() + nz; ++ i) {
        if(std::isnan(*i)) {
            *i = 0;
        }else if(std::isinf(*i) || std::isinf(-*i)) {
            *i = 0;
        }
    }
    Eigen::SparseMatrix<double> ret = - temp.transpose() * emPars.Mu[g];
    ret.makeCompressed();
    //cout << "Here" << endl;
    //cout << emPars.lambda << endl;
    //if(emPars.NE[g] > 2) {
    //    int st = emPars.NA[g] - emPars.NE[g] + 2;
    //    int len = emPars.NE[g] - 2;
    //    return *ret.valuePtr() + emPars.lambda * (emPars.A[g].block(st, 0, len, emPars.NX[g]) * X).sum();
    //}
    //return *ret.valuePtr() + emPars.lambda * X.block(emPars.NE[g], 0, emPars.NX[g] - emPars.NE[g], 1).sum();
    return *ret.valuePtr();
}

Eigen::MatrixXd Optimization::QGradient(Eigen::MatrixXd& X, int g) {
    //cout << "Now compute Q gradient" << endl;
    Eigen::SparseMatrix<double> sX = X.sparseView();
    Eigen::SparseMatrix<double> temp = emPars.C[g] * sX;
    temp.makeCompressed();
    temp = emPars.Mu[g].cwiseQuotient(temp);
    int nz = temp.nonZeros();
    for(auto i = temp.valuePtr(); i != temp.valuePtr() + nz; ++ i) {
        if(std::isnan(*i)) {
            *i = 0;
        }else if(std::isinf(*i) || std::isinf(-*i)) {
            *i = 0;
        }
    }
    
    //for(int w = 0 ; w < NW; w ++) {
    //    if(std::isnan(temp(w, 0))) {
    //        temp(w, 0) = 0;
    //    }
    //    if(std::isinf(temp(w, 0)) || std::isinf(-temp(w, 0))) {
    //        double s = 0;
    //        for(int i = 0; i < NX[g]; i ++) {
    //            s += C[g].coeffRef(w, i) * X(i, 0);
    //        }
    //        temp(w, 0) = 0;
    //    }
    //}
    temp = emPars.C[g].transpose() * temp;
    Eigen::MatrixXd jac = Eigen::MatrixXd(temp);
    //if(emPars.NE[g] >= 2) {
    //    int st = emPars.NA[g] - emPars.NE[g] + 2;
    //    int len = emPars.NE[g] - 2;
    //    auto aaaa = emPars.A[g].block(st, 0, len, emPars.NX[g]).colwise().sum().array() * emPars.lambda;
    //    jac = jac.array() + aaaa.transpose().array();
    //}
    //cout << "here" << endl; 
    //jac.block(emPars.NE[g], 0, emPars.NX[g] - emPars.NE[g], 1).array() += emPars.lambda;
    //jac.array() += emPars.lambda;

    //cout << "jac before " << jac << endl;
    //cout << "jac summation" << jac.sum() << endl;
    jac /= jac.norm();
    //cout << "jac " << jac << endl;
    ////cout << "jac size " << jac.rows() << ' ' << jac.cols() << endl;
    ////cout << "Summation is: " <<  jac.array().sum() << endl;
    Eigen::MatrixXd ret(emPars.NX[g] + 1, 1);
    ret.setZero(emPars.NX[g] + 1, 1);
    ret.block(0, 0, emPars.NX[g], 1) = -jac;
    return ret;
}

Eigen::MatrixXd Optimization::QConstraints(Eigen::MatrixXd& X, int g) {
    //cout << "Now compute Q constraints" << endl;
    Eigen::MatrixXd cons(emPars.NA[g] + 1, 1);
    //cons(0, 0) = X.sum() - 1;
    cons(0, 0) = (emPars.L[g] * X)(0, 0) - 1;
    cons.block(1, 0, emPars.NA[g], 1) = emPars.A[g] * X;
    //for(int i = 0 ; i < NA[g] + 1; i ++ ) {
    //    if(!(cons(i, 0) > -EPS)) {
    //        cout << "FALSE" << endl;
    //        cout << cons(i, 0) << endl;
    //    }
    //}
    return cons;
}

Eigen::MatrixXd Optimization::QConstraintsNormal(Eigen::MatrixXd& X, int g) {
    //cout << "Now compute Q constraints normal" << endl;
    Eigen::MatrixXd normals(emPars.NA[g] + 1, emPars.NX[g] + 1);
    normals.setZero(emPars.NA[g] + 1, emPars.NX[g] + 1);
    //normals.block(0, 0, 1, emPars.NX[g]) = Eigen::MatrixXd::Ones(1, emPars.NX[g]);
    normals.block(0, 0, 1, emPars.NX[g]) = emPars.L[g];
    normals.block(1, 0, emPars.NA[g], emPars.NX[g]) = emPars.A[g];
    //int cnt = 0;
    //for(int i = 0; i < NA[g] + 1; i ++) {
    //    for(int j = 0; j < NX[g] + 1; j ++) {
    //        normals(i, j) = cnt ++;
    //    }
    //}
    //return normals.transpose();
    return normals;
}

