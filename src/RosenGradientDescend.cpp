#include "RosenGradientDescend.h"

vector<double> RosenGradientDescend::timer(10, 0);

RosenGradientDescend::RosenGradientDescend(EMAlgorithm& emPars, int g): 
    emPars(emPars), NE(emPars.NE[g]), NA(emPars.NA[g]), NX(emPars.NX[g]), Njunc(NX - NE), g(g) { }

RosenGradientDescend::~RosenGradientDescend() { }

Eigen::MatrixXd RosenGradientDescend::optimizeQ(double eps, int iter) {
    this -> epsZero = 1e-12;
    //cout << endl << "Gene " << g << endl;
    //cout << "There are " << NE << " exons." << endl;
    //cout << "#Junctions = " << Njunc << endl;

    Eigen::MatrixXd X = emPars.X[g];
    Eigen::MatrixXd preX = X;

    //Eigen::MatrixXd ieq = Eigen::MatrixXd::Zero(NA + Njunc, NX);
    //ieq.block(0, 0, NA, NX) = emPars.A[g];
    //ieq.block(NA, NE, Njunc, Njunc) = Eigen::MatrixXd::Identity(Njunc, Njunc);
    //Eigen::SparseMatrix<double, Eigen::RowMajor> ieq(NA+Njunc, NX);
    //ieq.topRows(NA) = emPars.A[g];
    //for(int i = 0; i < Njunc; i ++) {
    //    ieq.insert(NA+i, NE+i) = 1.0;
    //}
    //ieq.makeCompressed();
    Eigen::SparseMatrix<double, Eigen::RowMajor> ieq(NA, NX);
    ieq.topRows(NA) = emPars.A[g];
    ieq.makeCompressed();
    //cout << ieq << endl;

    //Eigen::MatrixXd eq = Eigen::MatrixXd::Ones(1, NX);
    Eigen::MatrixXd eq = emPars.L[g];

    int it = 0;
    int smallStepCnt = 0;
    double t[10];
    memset(t, 0, sizeof(t));

    iter = (int)(NA * 1.2);
    Eigen::MatrixXd preD;

    while(it < iter){
        clock_t st = clock();
        //cout << endl;
        //cout << endl;
        //cout << "+++++++++++++++++++" << endl;
        //cout << "Gene " << g << " Iterations " << it << " includes #exons = " << NE << " #x = " << NX << endl;
        //cout << "QFunc = " << QFunction(X) << endl; 
        it ++;

        Eigen::MatrixXd cons = ieq * X;

//Eigen::internal::set_is_malloc_allowed(false);
        auto grad = QGradient(X);
        //if(grad.n)
        
//Eigen::internal::set_is_malloc_allowed(true);
        auto bnds = cons.array() < epsZero;

        //cout <<"X === " << endl;
        //cout << X.transpose() << endl;
        //cout <<"X ^^^ " << endl;

        //cout << "Cons === " << endl;
        //cout << cons.transpose() << endl;
        //cout << "Cons ^^^ " << endl;

        if(cons.minCoeff() < -epsZero) {
            //cout << "Gene " << g << " Iterations " << it << endl;
            cout << "Not Feasible!!" << cons.minCoeff() << endl;
            //getchar();
        }
        
        int Neqcons = 0, Nieqcons = 0;
        for(int i = 0; i < ieq.rows(); i ++) {
            bnds(i, 0) ? Neqcons ++ : Nieqcons ++;
        }

        //cout << "Number of hit " << Neqcons << endl;
        //cout << "Number of inner " << Nieqcons << endl;
        //Eigen::MatrixXd A1(Neqcons, NX);
        //Eigen::MatrixXd A2(Nieqcons, NX);
        Eigen::SparseMatrix<double, Eigen::RowMajor> A1(Neqcons, NX);
        Eigen::SparseMatrix<double, Eigen::RowMajor> A2(Nieqcons, NX);

        int keqcons = 0, kieqcons = 0;
        for(int i = 0; i < ieq.rows(); i ++) {
            if(bnds(i, 0)) {
                A1.middleRows<1>(keqcons ++) = ieq.row(i);
            } else {
                A2.middleRows<1>(kieqcons ++) = ieq.row(i);
            }
        }
        t[1] += (double)(clock() - st) / CLOCKS_PER_SEC;

        while(true) {
            //cout << endl << endl << endl;
            clock_t st2 = clock();
            //Eigen::MatrixXd M(A1.rows() + 1, NX);
            //M.middleRows(0, A1.rows()) = A1;
            //M.bottomRows<1>() = eq;

            //Eigen::SparseMatrix<double, Eigen::RowMajor> spM = M.sparseView(0.0);
            Eigen::SparseMatrix<double, Eigen::RowMajor> spM(A1.rows()+1, NX);
            spM.topRows(A1.rows()) = A1;
//Eigen::internal::set_is_malloc_allowed(false);
            spM.bottomRows<1>() = eq.sparseView();
//Eigen::internal::set_is_malloc_allowed(true);
            spM.makeCompressed();
            
            //cout << "X ===" << endl;
            //cout << X.sparseView() << endl;
            //cout << "X ^^^" << endl;

            //cout << "M sum = " << spM.sum() << " with " << spM.rows() << " rows" << endl;

            //cout << "Boundary CONS === " << endl;
            //for(int i = 0; i < ieq.rows(); i ++) {
            //    if(bnds(i, 0)) {
            //        cout << cons(i,0) << "\t";
            //    }
            //}
            //cout << endl;
            //cout << "Boundary CONS ^^^ " << endl;

            //cout << "M ===" << endl;
            //cout << spM << endl;
            //cout << "M ^^^" << endl;

            //cout << "Before inverse" << endl;
            clock_t st5 = clock();
            //Eigen::MatrixXd Inv = (M * M.transpose()).ldlt().solve(Eigen::MatrixXd::Identity(M.rows(), M.rows()));
            ////Eigen::MatrixXd Inv = (M * M.transpose()).inverse();
            //Eigen::MatrixXd R = Inv * M;

            //Eigen::FullPivLU <Eigen::MatrixXd> lu_decomp(spM);
            //auto rnk = lu_decomp.rank();
            //if(rnk < spM.rows()) {
            //    cout << "spM Rank " << rnk << " vs Rows " << spM.rows() << endl;
            //    //exit(-1);
            //}

            //Eigen::PardisoLLT<Eigen::SparseMatrix<double> > solver;
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
            //cout << "spM size " << spM.rows() << " " << spM.cols() << endl;
            //cout << "Inverse -- 1" << endl;
//Eigen::internal::set_is_malloc_allowed(false);
            Eigen::SparseMatrix<double> Q(spM * spM.transpose());
            solver.analyzePattern(Q);
            solver.factorize(Q);
//Eigen::internal::set_is_malloc_allowed(true);
            //solver.compute(spM * spM.transpose());
            //cout << "Inverse -- 2" << endl;
            Eigen::SparseMatrix<double> I(spM.rows(), spM.rows());
            //cout << "Inverse -- 3" << endl;
            I.setIdentity();
            //cout << "Inverse -- 4" << endl;
//Eigen::internal::set_is_malloc_allowed(false);
            Eigen::MatrixXd R = solver.solve(I) * spM;
//Eigen::internal::set_is_malloc_allowed(true);
            //R.array() += epsCtrl;
            //cout << "Inverse -- 5" << endl;

            //cout << "R sum = " << R.sum() << endl; 
            //cout << "After inverse" << endl;

            //cout << "||M*M'|| = " << (spM * spM.transpose()).norm() << endl;
            //cout << "M*M' ===" << endl;
            //cout << spM * spM.transpose() << endl;
            //cout << "M*M' ^^^" << endl;

            //cout << "M*M' -1 === " << endl;
            //cout << solver.solve(I) << endl;
            //cout << "M*M' -1 ^^^ " << endl;
            
            //cout << "R ===" << endl;
            //cout << R << endl;
            //cout << "R ^^^" << endl;

            t[5] += (double)(clock() - st5) / CLOCKS_PER_SEC;

            clock_t st6 = clock();
            //Eigen::MatrixXd P = Eigen::MatrixXd::Identity(NX, NX) - M.transpose() * R;

            //cout << "M' === " << endl;
            //cout << spM.transpose() << endl;
            //cout << "M' ^^^" << endl;


            //cout << "||R|| = " << R.norm() << endl;
            //cout << "M'*R === " << endl;
            //cout << spM.transpose() * R << endl;
            //cout << "M'*R ^^^ " << endl;

//Eigen::internal::set_is_malloc_allowed(false);
            Eigen::MatrixXd P = Eigen::MatrixXd::Identity(NX, NX) - spM.transpose() * R;
//Eigen::internal::set_is_malloc_allowed(true);

            //cout << "P sum = " << P.sum() << endl;
            t[6] += (double)(clock() - st6) / CLOCKS_PER_SEC;

            clock_t st7 = clock();
            //cout << "Checker" << " ||P|| = " << P.norm() << endl;
            //cout << "Checker" << " ||grad|| = " << grad.norm() << endl;
            

            Eigen::MatrixXd d = - P * grad;

            
            t[7] += (double)(clock() - st7) / CLOCKS_PER_SEC;
            t[2] += (double)(clock() - st2) / CLOCKS_PER_SEC;

            //cout << " P ===" << endl;
            //cout << P << endl;
            //cout << " P ^^^" << endl;

            //cout << "-grad-> ===== " << endl;
            //cout << grad.transpose() << endl;
            //cout << "-grad-> ^^^^^ " << endl;

            //cout << "-d-> ===== " << endl;
            //cout << d.transpose() << endl;
            //cout << "-d-> ^^^^^ " << endl;

            //cout << "||d|| = " << d.norm() << endl;
            //if(d.norm() - 1.0 > epsZero) {
            //    cout <<"NORMD over flow !!! " << endl;
            //    //exit(-1);
            //}
            

            if(d.norm() < epsZero) {
                clock_t st3 = clock();
                Eigen::MatrixXd W = R * grad;

                int mini = -1;
                double minv = -epsZero;

            /* NOTE HERE, here assume the eq must be the last row!!!!*/
                for(int i = 0; i < W.rows() - 1; i ++) {
                    if(W(i, 0) < minv) {
                        minv = W(i, 0);
                        mini = i;
                    }
                }
                t[3] += (double)(clock() - st3) / CLOCKS_PER_SEC;

                if(mini == -1) {
                    for(int i = 0; i < 8; i ++) {
                        cout << "T" << i << " = " << t[i] << endl;
                        timer[i] += t[i];
                    }
                    cout << "Total iteration is " << it << endl;
                    cout << "===Inner KKT============" << endl;
                    printf("%.6f\n", QFunction(X));
                    return X;
                } else {
                    int nRows = A1.rows() - 1;
                    int nCols = A1.cols();
                    A1.middleRows(mini, nRows-mini) = A1.middleRows(mini+1, nRows-mini);
                    A1.conservativeResize(nRows, nCols);
                    cout << "-------A row is removed-------" << endl;
                    continue;
                }
            } else {
                clock_t st4 = clock();

                d.array() /= d.norm();
                Eigen::MatrixXd dHat = A2 * d;
                Eigen::MatrixXd bHat = - A2 * X;
                //cout << "dCap === " << endl;
                //cout << dCap << endl;
                //cout << "dCap ^^^ " << endl;
                //
                //cout << "bCap === " << endl;
                //cout << bCap << endl;
                //cout << "bCap ^^^ " << endl;
                int dHatIdx = -1;
                double upperBound = 1.0;
                for(int i = 0; i < dHat.rows(); i ++) {
                    if(dHat(i, 0) < - epsZero) {
                        if(upperBound > bHat(i, 0) / dHat(i, 0)) {
                            upperBound = bHat(i, 0) / dHat(i, 0);
                            dHatIdx = i;
                        }
                        //cout << bCap(i, 0) << " / " << dCap(i, 0) << endl;
                        //cout << "Upper Bound === " << upperBound<< endl;
                    }
                }
                //cout << "UpperBound is " << upperBound << endl;
                //if(dHatIdx >= 0) {
                //    cout << "Corresponds to " << dHatIdx << " row and the cons is " << endl;
                //    cout << A2.row(dHatIdx) << endl; 
                //}
                double lambda = lineSearch(X, d, 0.0, upperBound);
                //lambda += epsCtrl;
                preX = X;

                X.noalias() += lambda * d;

                //cout << "Lambda is " << lambda << " and upper bound is " << upperBound << endl;
                //cout << "dX = " << (X-preX).norm() << endl;
                //cout << "QFunc " << QFunction(X) << endl;
                //cout << "X cons --> " << (ieq * X).minCoeff() << endl;
                //for(int i = 0; i < NX; i ++) {
                //    X(i, 0) = max(eps, min(1.0 - eps, X(i, 0)));
                //}
                //double tot = X.sum();
                //for(int i = 0; i < NX; i ++) {
                //    X(i, 0) /= tot;
                //}
                //cout << "QFunc " << QFunction(X) << endl;
                //cout << "X cons --> " << (ieq * X).minCoeff() << endl;

                t[4] += (double)(clock() - st4) / CLOCKS_PER_SEC;
                if(lambda < epsZero) {
                    for(int i = 0; i < 8; i ++) {
                        cout << "T" << i << " = " << t[i] << endl;
                        timer[i] += t[i];
                    }
                    cout << "Total iteration is " << it << endl;
                    cout << "===Small Step Length============" << endl;
                    printf("%.6f\n", QFunction(X));
                    return X;
                }/* else {
                    if(lambda < epsZero) {
                        smallStepCnt ++;
                        if(smallStepCnt > 20) {
                            cout << "Total iteration is " << it << endl;
                            cout << "===Much Small Step Length============" << endl;
                            return X;
                        }
                    }
                    break;
                }*/
                break;
            }
        }
        
        t[0] += double(clock() - st) / CLOCKS_PER_SEC;
    }
    for(int i = 0; i < 8; i ++) {
        cout << "T" << i << " = " << t[i] << endl;
        timer[i] += t[i];
    }
    cout << "Total iteration is " << it << endl;
    cout << "===Iteration exceed ============" << endl;
    printf("%.6f\n", QFunction(X));
    return X;
}

double RosenGradientDescend::lineSearch(Eigen::MatrixXd& X, Eigen::MatrixXd& d, double lowerBound, double upperBound) {
    //cout << "Line search X === " << endl;
    //cout << X << endl;
    //cout << "Line search X ^^^ " << endl;
    double l = lowerBound, r = upperBound;
    double lm = l + (1.0 - GoldenSplit) * (r - l), rm = l + GoldenSplit * (r - l);
    Eigen::MatrixXd Xlm = X;
    Xlm.noalias() += lm * d;
    //cout << "Line search XLM === " << endl;
    //cout << Xlm << endl;
    //cout << "Line search XLM ^^^ " << endl;
    Eigen::MatrixXd Xrm = X;
    Xrm.noalias() += rm * d;
    //cout << "Line search XRM === " << endl;
    //cout << Xrm << endl;
    //cout << "Line search XRM ^^^ " << endl;
    //cout << "LM = " << lm << endl;
    double lf = QFunction(Xlm);
    //cout << "RM = " << rm << endl;
    double rf = QFunction(Xrm);
    double epsF = pow(10, (int)(log(lf) / log(10.0)) - 13);
    //cout << "epsF = " << epsF << endl;
    
    while(l - r < - epsZero) {
        //cout << l - r << endl;
        //cout << "Line Search: " << lf + epsZero << " vs " << rf + epsZero << "\t" << l << " vs " << r << endl;
        //cout << rf - lf << endl;
        if(lf - rf < - epsF) {
            r = rm - epsZero;
            rm = lm;
            rf = lf;
            lm = l + (1 - GoldenSplit) * (r - l);
            Xlm = X;
            Xlm.noalias() += lm * d;
            lf = QFunction(Xlm);
        } else {
            l = lm + epsZero;
            lm = rm;
            lf = rf;
            rm = l + GoldenSplit * (r - l);
            Xrm = X;
            Xrm.noalias() += rm * d;
            rf = QFunction(Xrm);
        }
    }
    return r;
}

double RosenGradientDescend::QFunction(Eigen::MatrixXd& X) {
    Eigen::SparseMatrix<double> sX = X.sparseView();
    //cout << "Inner dense X ===" << endl; 
    //cout << X << endl;
    //cout << "Inner dense X ^^^" << endl; 
    sX.makeCompressed();
    //cout << "Inner sparse X === " << endl;
    //cout << X << endl;
    //cout << "Inner sparse X ^^^ " << endl;
    Eigen::SparseMatrix<double> temp = emPars.C[g] * sX;
    temp.makeCompressed();
    //cout << "||Inner temp before log|| = " << temp.norm() << endl;
    //cout << "Inner temp before log === " << endl;
    //cout << temp << endl;
    //cout << "Inner temp before log ^^^ " << endl;
    int nz = temp.nonZeros();
    for(auto i = temp.valuePtr(); i != temp.valuePtr() + nz; ++ i) {
        *i = log(*i);
    }
    //cout << "||temp|| = " << temp.norm() << endl;
    //cout << "Inner temp after log === " << endl;
    //cout << temp << endl;
    //cout << "Inner temp after log ^^^ " << endl;
    Eigen::SparseMatrix<double> ret = - temp.transpose() * emPars.Mu[g];
    return *ret.valuePtr();
}

Eigen::MatrixXd RosenGradientDescend::QGradient(Eigen::MatrixXd& X) {
    Eigen::SparseMatrix<double> sX = X.sparseView();
    sX.makeCompressed();
    Eigen::SparseMatrix<double> temp = emPars.C[g] * sX;
    temp.makeCompressed();
    //cout << "Compressed? " << temp.isCompressed() << endl;
    //cout << "Checktemp-1 " << temp.norm() << endl;
    //cout << temp.nonZeros() << endl;
    //cout << "TEMPTEMP" << endl;
    //cout << temp << endl;
    //cout << "MUMU" << endl;
    //cout << emPars.Mu[g] << endl;
    //EigenSparseMatrix<double> tempMask = temp.cwiseQuotient(temp);
    temp = emPars.Mu[g].cwiseProduct(temp.cwiseInverse());
    temp.makeCompressed();

    //for(auto i = temp.valuePtr(); i != temp.valuePtr() + temp.nonZeros(); ++ i) {
    //    if(std::isnan(*i) || std::isinf(*i)) {
    //        *i = 0;
    //    }
    //}

    //cout << "Checktemp-2 " << temp.norm() << endl;
    //cout << temp.nonZeros() << endl;
    //cout << temp << endl;
    //cout << "CheckMu-2 " << emPars.Mu[g].norm() << endl;
//Eigen::internal::set_is_malloc_allowed(false);
    temp = emPars.C[g].transpose() * temp;
//Eigen::internal::set_is_malloc_allowed(true);

    //cout << "Checktemp-3 " << temp.norm() << endl;
    temp /= - temp.norm();
    return temp;
    //Eigen::MatrixXd jac = temp;
    //if(std::isnan(jac.norm())) {
    //    cout << "!!!!FATAL jac!!!!" << endl;
    //    cout << X << endl;
    //    cout << "!!!!FATAL jac!!!!" << endl;
    //    cout << jac << endl;
    //    cout << "!!!!FATAL SX!!!!" << endl;
    //    cout << sX << endl;
    //    cout << "!!!!FATAL temp!!!!" << endl;
    //    cout << temp << endl;
    //}
    //cout << jac.norm() << " inner grad" << endl;
    //cout << jac.sum() << " sum inner grad" << endl;
    //jac /= jac.norm();
    //return - jac;
}

