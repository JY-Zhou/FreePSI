#include "ConjugateGradientProjection.h"

vector<double> ConjugateGradientProjection::timer(25, 0);
vector<double> ConjugateGradientProjection::exitPort(10, 0);
double ConjugateGradientProjection::QFunc_t = 0.0;
double ConjugateGradientProjection::QGrad_t = 0.0;
double ConjugateGradientProjection::ConsAdd_t = 0.0;
int ConjugateGradientProjection::QFunc_n = 0;
int ConjugateGradientProjection::QGrad_n = 0;
int ConjugateGradientProjection::ConsAdd_n = 0;
int ConjugateGradientProjection::iteration = 0;

ConjugateGradientProjection::ConjugateGradientProjection(EMAlgorithm& emPars, int g): 
    emPars(emPars), NX(emPars.NX[g]), NEQ(1), g(g), epsZero(1e-9) { }

ConjugateGradientProjection::~ConjugateGradientProjection() { }

void ConjugateGradientProjection::initOptimization() {
    /*** Initialze constraints ***/
    ieq = emPars.A[g].transpose();
    //eq = emPars.L[g].transpose();
    eq = Eigen::MatrixXd::Ones(NX, 1);
    for(int j = 0; j < ieq.cols(); j ++) {
        ieq.col(j) /= ieq.col(j).norm();
    }
    eq.array() /= eq.norm();

    /*** Initialze variables ***/
    CC = emPars.C[g].transpose() * emPars.C[g];
    CM = emPars.C[g].transpose() * emPars.Mu[g];
    X = Eigen::MatrixXd(emPars.X[g]);
    grad = QGradient(X);
    
    /*** Initialze constraint indices ***/
    hyperplaneSeq.clear();
    stateOfIeqs.clear();
    for(int j = 0; j < ieq.cols(); j ++) {
        stateOfIeqs.push_back(-1);
    }

    /*** Initialze INV, P and H ***/
    hyperplaneSeq.push_back(-1);
    INV = Eigen::MatrixXd::Identity(1, 1);
    N = eq.sparseView();
    //P = Eigen::MatrixXd::Identity(NX, NX) - eq * eq.transpose();
    P = Eigen::MatrixXd::Identity(NX, NX);
    for(int j = 0; j < P.cols(); j ++) {
        P.col(j).noalias() -= eq(j, 0) * eq;
    }
    //H.resize(NX, NX);
    //H.setIdentity();
    //double temp = (eq.transpose() * H * eq)(0, 0);
    //H = H - (H * N * N.transpose() * H) / temp;
    H = Eigen::MatrixXd::Identity(NX, NX);
    Eigen::MatrixXd Heq = H * eq;
    double eqHeq = Heq.col(0).dot(eq.col(0));
    for(int j = 0; j < H.cols(); j ++) {
        H.col(j).noalias() -= Heq(j, 0) / eqHeq * Heq;
    }

    /*** Check the activity of all constrants ***/
    for(int j = 0; j < ieq.cols(); j ++) {
        //if(Eigen::SparseMatrix<double>(ieq.col(j).transpose() * X).coeffRef(0, 0) < epsZero) {
        if(ieq.col(j).dot(X.col(0)) < epsZero) {
            if(addHyperplane(j)) {
                stateOfIeqs[j] = hyperplaneSeq.size() - NEQ;
                hyperplaneSeq.push_back(j);
            }
        }
    }
}

Eigen::MatrixXd ConjugateGradientProjection::optimizeQ(double eps, int iter) {
    clock_t st0 = clock();
    initOptimization();
    timer[1] += (double)(clock() - st0) / CLOCKS_PER_SEC;
    //iteration = 0;

    for(int it = 0; it < iter; it ++) {
        /*** Step 1 ***/
        clock_t st1 = clock();
        //cout << "step 1" << endl;
        Eigen::MatrixXd s = H * grad;
        //cout << "H === " << H << endl;
        //cout << "grad === " << grad << endl;
        //cout << "s === " << s << endl;
        Eigen::MatrixXd alpha = INV * (N.transpose() * grad);
        bool checkFlag = true;
        for(int i = NEQ; i < alpha.rows(); i ++) {
            if(alpha(i, 0) > epsZero) {
                checkFlag = false;
            } 
        }

        double sLen = s.norm();
        if(sLen < 1e-8 && checkFlag) {
            //cerr << "*" << g << "**** Constrained Stationary Point *****" << endl; 
            //cout << "***** Constrained Stationary Point *****" << endl; 
            timer[0] += (double)(clock() - st0) / CLOCKS_PER_SEC;
            //cout << timer[0] << endl;
            exitPort[0] += 1;
            iteration += it + 1;
            return X;
        }
        timer[2] += (double)(clock() - st1) / CLOCKS_PER_SEC;

        /*** Step 2 ***/
        clock_t st2 = clock();
        //cout << "step 2" << endl;
        double maxAlphaB = -epsZero;
        int maxAlphaBIdx = -1;
        for(int i = NEQ; i < alpha.rows(); i ++) {
            double alphaB = 0.5 * alpha(i, 0) / sqrt(INV(i, i));
            if(maxAlphaB - alphaB < epsZero) {
                maxAlphaB = alphaB;
                maxAlphaBIdx = i;
            }
        }
        timer[3] += (double)(clock() - st2) / CLOCKS_PER_SEC;

        if(sLen - maxAlphaB < epsZero) {
            /*** drop the last hyperplane ***/
            clock_t st3 = clock();
            //Eigen::SparseMatrix<double> q = N.col(maxAlphaBIdx);
            swapLastHyperplane(maxAlphaBIdx);
            dropLastHyperplane();
            stateOfIeqs[hyperplaneSeq.back()] = -1;
            hyperplaneSeq.pop_back();
            timer[4] += (double)(clock() - st3) / CLOCKS_PER_SEC;
        } else {
            /*** Step 3 ***/
            clock_t st4 = clock();
            //cout << "step 3" << endl;
            s /= sLen;

            double upperBound = 1e3;
            int upperBoundIdx = -1;
            vector <double> zeroBoundIdx; 
            for(int j = 0; j < ieq.cols(); j ++) {
                if(stateOfIeqs[j] < 0) {
                    //double sProj = Eigen::SparseMatrix<double>(ieq.col(j).transpose() * s).coeffRef(0, 0);
                    double sProj = ieq.col(j).dot(s.col(0));
                    if(fabs(sProj) < epsZero) {
                        //cout << " is the linear dependency one. PASS!" << endl;
                    } else {
                        //double bound = - Eigen::SparseMatrix<double>(ieq.col(j).transpose() * X).coeffRef(0, 0) / sProj;
                        double bound = - ieq.col(j).dot(X.col(0)) / sProj;

                        if(fabs(bound) < epsZero) {
                            zeroBoundIdx.push_back(j);
                        } else if (bound > epsZero) {
                            if(bound - upperBound < epsZero) {
                                upperBound = bound;
                                upperBoundIdx = j;
                            }
                        }
                    }
                }
            }
            timer[5] += (double)(clock() - st4) / CLOCKS_PER_SEC;

            /*** Do line search ***/
            clock_t st5 = clock();
            //double lambda = GoldenSelection(X, s, 0.0, upperBound);
            double lambda = SecantMethod(X, s, 0.0, upperBound);
            //cout << g << " === " << lambda << endl; 
            //cout << s << endl;
            //cout << "~~" << endl;
            //cout << X << endl;
            //getchar();

            if(lambda - 1e3 > -epsZero) {
                exitPort[1] += 1;
                timer[0] += (double)(clock() - st0) / CLOCKS_PER_SEC;
                iteration += it + 1;
                return X;
            }

            //Eigen::SparseMatrix<double> newX = X + lambda * s;
            Eigen::MatrixXd newX = X + lambda * s;
            //cout << "XNZ = " <<  newX.nonZeros() << endl;
            //Eigen::SparseMatrix<double> newGrad = QGradient(newX);
            Eigen::MatrixXd newGrad = QGradient(newX);
            //cout << "Grad NZ = " <<  newX.nonZeros() << endl;
            timer[6] += (double)(clock() - st5) / CLOCKS_PER_SEC;

            /*** Check new constraints ***/
            clock_t st6 = clock();
            //Eigen::SparseMatrix<double> checkIEQ = ieq.transpose() * newX;
            bool addNewBound = false;
            for(int j = 0; j < ieq.cols(); j ++) {
                //if(checkIEQ.coeffRef(j, 0) < - epsZero) {
                if(ieq.col(j).dot(newX.col(0)) < -epsZero) {
                    if(find(zeroBoundIdx.begin(), zeroBoundIdx.end(), j) != zeroBoundIdx.end()) {
                        if(addHyperplane(j)) {
                            stateOfIeqs[j] = hyperplaneSeq.size() - NEQ;
                            hyperplaneSeq.push_back(j);
                            addNewBound = true;
                        }
                    } else {
                        //newX = correctInfeasible(newX);
                    }
                }
            }

            if(lambda < epsZero) {
                //cout << "***** Small Step Length *****" << endl;
                exitPort[2] += 1;
                timer[0] += (double)(clock() - st0) / CLOCKS_PER_SEC;
                iteration += it + 1;
                return X;
            }
            timer[7] += (double)(clock() - st6) / CLOCKS_PER_SEC;

            if(!addNewBound) {
                /*** Step 4 ***/
                //cout << "step 4" << endl;
                if(fabs(lambda - upperBound) < epsZero) {
                    //cout << "***** Lambda touch the bound *****" << endl;
                    clock_t st7 = clock();
                    if(addHyperplane(upperBoundIdx)) {
                        stateOfIeqs[upperBoundIdx] = hyperplaneSeq.size() - NEQ;
                        hyperplaneSeq.push_back(upperBoundIdx);
                    }
                    timer[8] += (double)(clock() - st7) / CLOCKS_PER_SEC;
                } else {
                    //cout << "***** Update H because X is going to be updated *****" << endl;
                    clock_t st8 = clock();
                    //Eigen::SparseMatrix<double> sigma = lambda * s;
                    Eigen::MatrixXd sigma = lambda * s;
                    //Eigen::SparseMatrix<double> y = newGrad - grad;
                    Eigen::MatrixXd y = newGrad - grad;
                    //Update H
                    //Eigen::SparseMatrix<double> A = - sigma * sigma.transpose();
                    //A /= Eigen::SparseMatrix<double>(sigma.transpose() * y).coeffRef(0, 0);
                    //Eigen::SparseMatrix<double> B = - (H * y) * (y.transpose() * H);
                    //B /= Eigen::SparseMatrix<double>(y.transpose() * H * y).coeffRef(0, 0);
                    //H = H + A + B;
                    Eigen::MatrixXd Hy = H * y;
                    double yHy = Hy.col(0).dot(y.col(0));
                    double sy = sigma.col(0).dot(y.col(0));
                    for(int j = 0; j < H.cols(); j ++) {
                        H.col(j).noalias() -= sigma(j, 0) / sy * sigma + Hy(j, 0) / yHy * Hy;
                    }
                    
                    timer[9] += (double)(clock() - st8) / CLOCKS_PER_SEC;
                }

                X = newX;
                grad = newGrad;
            }
        }
    }
    //cout << "***** Iteration overflow *****" << endl;
    exitPort[3] += 1;
    timer[0] += (double)(clock() - st0) / CLOCKS_PER_SEC;
    iteration += iter;
    return X;
}

void ConjugateGradientProjection::swapLastHyperplane(int idxOfN) {
    //cout << "swap" << endl;
    if(idxOfN == N.cols() - 1) {
        return;
    }

    //Swap N
    Eigen::SparseMatrix<double> tempN = N.col(idxOfN);
    N.col(idxOfN) = N.rightCols<1>();
    N.rightCols<1>() = tempN;

    //Swap INV
    Eigen::MatrixXd tempINV = INV.row(idxOfN);
    INV.row(idxOfN) = INV.bottomRows<1>();
    INV.bottomRows<1>() = tempINV;
    tempINV = INV.col(idxOfN);
    INV.col(idxOfN) = INV.rightCols<1>();
    INV.rightCols<1>() = tempINV;

    //Swap hyperplaneSeq
    int i = hyperplaneSeq[idxOfN];
    int q = hyperplaneSeq[N.cols() - 1];
    hyperplaneSeq[idxOfN] = q;
    hyperplaneSeq[N.cols() - 1] = i;
    //cout << i << " swap with " << q << endl;

    //Swap stateOfIeqs
    int tempState = stateOfIeqs[i];
    stateOfIeqs[i] = stateOfIeqs[q];
    stateOfIeqs[q] = tempState;
}

void ConjugateGradientProjection::dropLastHyperplane() {
    //cout << "***** Remove the last hyperplane from N *****" << endl;
    //Compute (N_{q-1}'N_{q-1})^{-1}
    int nR = INV.rows();
    int nC = INV.cols();
    Eigen::MatrixXd B2 = INV.topRightCorner(nR - 1, 1);
    double B4 = 1.0 / INV(nR - 1, nC - 1);
    INV.conservativeResize(nR - 1, nC - 1);
    //INV = INV - B2 * B4 * B2.transpose();
    for(int j = 0; j < INV.cols(); j ++) {
        INV.col(j).noalias() -= B4 * B2(j, 0) * B2;
    }

    //Compute N_{q-1}
    Eigen::SparseMatrix<double> q = N.rightCols<1>();
    N.conservativeResize(NX, N.cols() - 1);

    //Compute P_{q-1} :: Maybe overhead!
    P = Eigen::MatrixXd::Identity(NX, NX) - N * INV * N.transpose();

    //Compute H_q
    //double temp = (q.transpose() * P * q)(0, 0);
    //Eigen::SparseMatrix<double> sparP = P.sparseView();
    //H = H + sparP * q * q.transpose() * sparP / temp;
    Eigen::MatrixXd Pq = P * q;
    double qPq = q.col(0).dot(Pq.col(0));
    for(int j = 0; j < H.cols(); j ++) {
        H.col(j).noalias() += Pq(j, 0) / qPq * Pq;
    }
}

bool ConjugateGradientProjection::addHyperplane(int j) {
    clock_t st = clock();
    clock_t st1 = clock();
    //cout << "***** Add Hyperplane IEQ_idx = " << j << " *****" << endl;
    Eigen::SparseMatrix<double> q = ieq.col(j);
    Eigen::MatrixXd r = INV * (N.transpose() * q);
    Eigen::MatrixXd u = Eigen::MatrixXd(q) - N * r;
    if(u.norm() < epsZero) {
        //cout << "=== Linear Independent Violated ===" << endl;
        //exit(-1);
        return false;
    }
    timer[10] += (double)(clock() - st1)/CLOCKS_PER_SEC;

    clock_t st2 = clock();
    //Compute (N_q'N_q)^{-1}
    double A0 = u.norm();
    A0 = 1.0 / (A0 * A0);
    //Eigen::MatrixXd B1 = INV + A0 * r * r.transpose();
    Eigen::MatrixXd B1 = INV;
    for(int j = 0; j < B1.cols(); j ++) {
        B1.col(j).noalias() +=  A0 * r(j, 0) * r;
    }
    Eigen::MatrixXd B2 = - A0 * r;
    int nR = INV.rows();
    int nC = INV.cols();
    INV.conservativeResize(nR + 1, nC + 1);
    INV.topLeftCorner(nR, nC) = B1;
    INV.topRightCorner(nR, 1) = B2;
    INV.bottomLeftCorner(1, nC) = B2.transpose();
    INV(nR, nC) = A0;
    timer[11] += (double)(clock() - st2)/CLOCKS_PER_SEC;

    clock_t st3 = clock();
    //Compute P_q
    //P = P - A0 * (u * u.transpose());
    for(int j = 0; j < P.cols(); j ++) {
        P.col(j).noalias() -= A0 * u(j, 0) * u;
    }
    timer[12] += (double)(clock() - st3)/CLOCKS_PER_SEC;

    clock_t st4 = clock();
    //Compute N_q
    N.conservativeResize(NX, N.cols() + 1);
    N.rightCols<1>() = q;
    timer[13] += (double)(clock() - st4)/CLOCKS_PER_SEC;

    clock_t st5 = clock();
    //Compute H_q
    //double temp = Eigen::SparseMatrix<double>(q.transpose() * H * q).coeffRef(0, 0);
    timer[14] += (double)(clock() - st5)/CLOCKS_PER_SEC;

    clock_t st6 = clock();
    //H = H - (H * q) * (q.transpose() * H) / temp;
    Eigen::MatrixXd Hq = H * q;
    double qHq = q.col(0).dot(Hq.col(0));
    timer[15] += (double)(clock() - st6)/CLOCKS_PER_SEC;

    clock_t st7 = clock();
    //Eigen::SparseMatrix<double> HqqH = Hq * Hq.transpose();
    timer[16] += (double)(clock() - st7)/CLOCKS_PER_SEC;

    clock_t st8 = clock();
    //H -= HqqH / temp;
    //for(Eigen::SparseMatrix<double>::InnerIterator it(Hq, 0); it; ++ it) {
    //    H.col(it.row()) -= Hq * Hq.coeffRef(it.row(), it.col()) / temp;
    //}
    for(int j = 0; j < H.cols(); j ++) {
        H.col(j).noalias() -= Hq(j, 0) / qHq * Hq;
    }
    //cout << H.rows() << " " << H.cols() << " " << H.nonZeros() << endl;
    timer[17] += (double)(clock() - st8)/CLOCKS_PER_SEC;
    ConsAdd_t += (double)(clock() - st) / CLOCKS_PER_SEC;
    ConsAdd_n ++;
    return true;
}

Eigen::SparseMatrix<double> ConjugateGradientProjection::correctInfeasible(Eigen::SparseMatrix<double> & X) {
    cout << "Poor!" << endl;
    //Project to equalities
    Eigen::SparseMatrix<double> Nvio = eq.sparseView();
    Eigen::MatrixXd INVvio = Eigen::MatrixXd::Identity(1, 1);
    Eigen::SparseMatrix<double> b;

    //cout << X << endl;
    //cout << "^^^ X" << endl;
    b.resize(1, 1), b.insert(0, 0) = 1;
    X = X - Nvio * (Nvio.transpose() * X - b);

    //Project to inequalities
    Eigen::SparseMatrix<double> checkIEQ = ieq.transpose() * X;
    vector <int> unchecked;
    for(int j = 0; j < ieq.cols(); j ++) {
        if(checkIEQ.coeffRef(j, 0) < epsZero) {
            Eigen::SparseMatrix<double> q = ieq.col(j);
            Eigen::MatrixXd r = INVvio * (Nvio.transpose() * q);
            Eigen::MatrixXd u = Eigen::MatrixXd(q) - Nvio * r;
            if(u.norm() < epsZero) {
                continue;
            }

            //Compute (N_q'N_q)^{-1}
            double A0 = u.norm();
            A0 = 1.0 / (A0 * A0);
            Eigen::MatrixXd B1 = INVvio + A0 * r * r.transpose();
            Eigen::MatrixXd B2 = - A0 * r;
            int nR = INVvio.rows();
            int nC = INVvio.cols();
            INVvio.conservativeResize(nR + 1, nC + 1);
            INVvio.topLeftCorner(nR, nC) = B1;
            INVvio.topRightCorner(nR, 1) = B2;
            INVvio.bottomLeftCorner(1, nC) = B2.transpose();
            INVvio(nR, nC) = A0;

            //Compute N_q
            Nvio.conservativeResize(NX, Nvio.cols() + 1);
            Nvio.rightCols<1>() = q;
        } else {
            unchecked.push_back(j);
        }
    }

    b.resize(Nvio.cols(), 1), b.insert(0, 0) = 1;
    Eigen::SparseMatrix<double> sparINVvio = INVvio.sparseView();
    X = X - Nvio * sparINVvio * (Nvio.transpose() * X - b);


    int pres = unchecked.size() + 1;

    while(!unchecked.empty()) {
        if(pres == unchecked.size()) {
            cout <<"Error here !! no hyperplane removed" << endl;
            exit(-1);
        }
        pres = unchecked.size();
        double minv = 1e9;
        int mini = -1;
        for(int i = 0; i < unchecked.size(); i ++) {
            int j = unchecked[i];
            Eigen::SparseMatrix<double> res = ieq.col(j).transpose() * X;
            if(minv - res.coeffRef(0, 0) > epsZero) {
                minv = res.coeffRef(0, 0);
                mini = i;
            }
        }

        if(minv > - epsZero) {
            return X;
        }

        //Projected to minimum
        Eigen::SparseMatrix<double> q = ieq.col(unchecked[mini]);
        unchecked.erase(unchecked.begin() + mini);

        Eigen::MatrixXd r = INVvio * (Nvio.transpose() * q);
        Eigen::MatrixXd u = Eigen::MatrixXd(q) - Nvio * r;

        if(u.norm() > epsZero) {
            //cout << "State 2" << endl;
            //Compute (N_q'N_q)^{-1}
            double A0 = u.norm();
            A0 = 1.0 / (A0 * A0);
            Eigen::MatrixXd B1 = INVvio + A0 * r * r.transpose();
            Eigen::MatrixXd B2 = - A0 * r;
            int nR = INVvio.rows();
            int nC = INVvio.cols();
            INVvio.conservativeResize(nR + 1, nC + 1);
            INVvio.topLeftCorner(nR, nC) = B1;
            INVvio.topRightCorner(nR, 1) = B2;
            INVvio.bottomLeftCorner(1, nC) = B2.transpose();
            INVvio(nR, nC) = A0;
            //Compute N_q
            Nvio.conservativeResize(NX, Nvio.cols() + 1);
            Nvio.rightCols<1>() = q;

            Eigen::SparseMatrix<double> res(Nvio.cols(), 1);
            res.insert(res.rows()-1, 0) = Eigen::SparseMatrix<double>(q.transpose() * X).coeffRef(0, 0);
            Eigen::SparseMatrix<double> sparinvvio = INVvio.sparseView();

            X = X - Nvio * sparinvvio * res;
        } else {
            //cout << "State 3" << endl;
            bool allMinus = true;
            for(int i = 0; i < r.rows(); i ++) {
                if(r(i, 0) > epsZero) {
                    allMinus = false;

                    //Swap with last one
                    if(i != Nvio.cols() - 1) {
                        //Swap N
                        Eigen::SparseMatrix<double> tempN = Nvio.col(i);
                        Nvio.col(i) = Nvio.rightCols<1>();
                        Nvio.rightCols<1>() = tempN;
                        //Swap INV
                        Eigen::MatrixXd tempINV = INVvio.row(i);
                        INVvio.row(i) = INVvio.bottomRows<1>();
                        INVvio.bottomRows<1>() = tempINV;
                        tempINV = INVvio.col(i);
                        INVvio.col(i) = INVvio.rightCols<1>();
                        INVvio.rightCols<1>() = tempINV;
                    }

                    if(true) {
                        //Drop last one
                        //Compute (N_{q-1}'N_{q-1})^{-1}
                        int nR = INVvio.rows();
                        int nC = INVvio.cols();
                        Eigen::MatrixXd B2 = INVvio.topRightCorner(nR - 1, 1);
                        double B4 = 1.0 / INVvio(nR - 1, nC - 1);
                        INVvio.conservativeResize(nR - 1, nC - 1);
                        INVvio = INVvio - B2 * B4 * B2.transpose();
                        //Compute N_{q-1}
                        Nvio.conservativeResize(NX, Nvio.cols() - 1);
                    }
                    r = INVvio * (Nvio.transpose() * q);
                    u = Eigen::MatrixXd(q) - Nvio * r;

                    if(true) {
                        //Add new
                        //Compute (N_q'N_q)^{-1}
                        double A0 = u.norm();
                        A0 = 1.0 / (A0 * A0);
                        Eigen::MatrixXd B1 = INVvio + A0 * r * r.transpose();
                        Eigen::MatrixXd B2 = - A0 * r;
                        int nR = INVvio.rows();
                        int nC = INVvio.cols();
                        INVvio.conservativeResize(nR + 1, nC + 1);
                        INVvio.topLeftCorner(nR, nC) = B1;
                        INVvio.topRightCorner(nR, 1) = B2;
                        INVvio.bottomLeftCorner(1, nC) = B2.transpose();
                        INVvio(nR, nC) = A0;
                        //Compute N_q
                        Nvio.conservativeResize(NX, Nvio.cols() + 1);
                        Nvio.rightCols<1>() = q;
                    }

                    Eigen::SparseMatrix<double> res(Nvio.cols(), 1);
                    res.insert(res.rows()-1, 0) = Eigen::SparseMatrix<double>(q.transpose() * X).coeffRef(0, 0);
                    Eigen::SparseMatrix<double> sparinvvio = INVvio.sparseView();

                    X = X - Nvio * sparinvvio * res;
                }
            }

            if(allMinus) {
                cout << "Error: Feasible space wont be empty!" << endl;
                exit(-1);
            }
        }
    }

    return X;
}

double ConjugateGradientProjection::GoldenSelection(Eigen::SparseMatrix<double>& X, Eigen::SparseMatrix<double>& d, double lowerBound, double upperBound) {
    double l = lowerBound, r = upperBound;
    double lm = l + (1.0 - GoldenSplit) * (r - l), rm = l + GoldenSplit * (r - l);
    Eigen::SparseMatrix<double> Xlm = X + lm * d;
    Eigen::SparseMatrix<double> Xrm = X + rm * d;
    double lf = QFunction(Xlm);
    double rf = QFunction(Xrm);
    double epsF = 1e-8;

    while(l - r < - epsZero) {
        if(rf - lf < epsF) {
            r = rm - epsZero;
            rm = lm;
            rf = lf;
            lm = l + (1 - GoldenSplit) * (r - l);
            Xlm = X + lm * d;
            lf = QFunction(Xlm);
        } else {
            l = lm + epsZero;
            lm = rm;
            lf = rf;
            rm = l + GoldenSplit * (r - l);
            Xrm = X + rm * d;
            rf = QFunction(Xrm);
        }
    }
    return r;
}

double ConjugateGradientProjection::SecantMethod(Eigen::MatrixXd & X, Eigen::MatrixXd & d, double lowerBound, double upperBound) {
    double a = lowerBound + (1.0 - GoldenSplit) * (upperBound - lowerBound);
    double b = lowerBound + GoldenSplit * (upperBound - lowerBound);
    int it = 0;
    while(fabs(a - b) > 1e-8 && it < 30) {
        //Eigen::SparseMatrix<double> Xa = X + a * d, Xb = X + b * d;
        Eigen::MatrixXd Xa = X + a * d, Xb = X + b * d;
        //double da = Eigen::SparseMatrix<double>(QGradient(Xa).transpose() * d).coeffRef(0, 0);
        //double db = Eigen::SparseMatrix<double>(QGradient(Xb).transpose() * d).coeffRef(0, 0);
        double da = d.col(0).dot(QGradient(Xa).col(0));
        double db = d.col(0).dot(QGradient(Xb).col(0));
        double c = b - (b-a)/(db-da)*db;
        c = min(upperBound, max(lowerBound, c));
        if(c == a) {
            return max(a, b);
        }
        a = b;
        b = c;
        it ++;
    }
    return b;
}

double ConjugateGradientProjection::QFunction(Eigen::SparseMatrix<double>& X) {
    clock_t st = clock();
    X.makeCompressed();
    Eigen::SparseMatrix<double> temp = emPars.C[g] * X;
    //cout << endl;
    //cout << "CX norm = " << temp.norm() << endl;
    temp.prune(epsZero, 0.01);
    temp.makeCompressed();
    int nz = temp.nonZeros();
    //cout << nz << endl;
    for(auto i = temp.valuePtr(); i != temp.valuePtr() + nz; ++ i) {
        //if(*i < epsZero)
            //cout << *i << " ";
        *i = log(*i);
    }
    //cout << endl;
    //cout << "CX norm = " << temp.norm() << endl;
    //cout << Eigen::SparseMatrix<double>(emPars.Mu[g].transpose() * temp).coeffRef(0, 0) << endl;
    double ans = Eigen::SparseMatrix<double>(emPars.Mu[g].transpose() * temp).coeffRef(0, 0);
    QFunc_t += (double)(clock() - st) / CLOCKS_PER_SEC;
    QFunc_n ++;
    return ans;
}

Eigen::MatrixXd ConjugateGradientProjection::QGradient(Eigen::MatrixXd & X) {
    clock_t st = clock();

    clock_t st1 = clock();
    Eigen::SparseMatrix<double> sX = X.sparseView();
    timer[18] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    //st1 = clock();
    //cout << emPars.C[g].rows() << " " << emPars.C[g].cols() << " " << emPars.C[g].nonZeros() << endl;
    //Eigen::SparseMatrix<double> CX = emPars.C[g] * sX;
    //timer[19] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    //st1 = clock();
    //CX.prune(epsZero, 0.01);
    //CX.makeCompressed();
    //timer[20] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    //st1 = clock();
    //Eigen::MatrixXd temp = emPars.Mu[g].cwiseProduct(CX.cwiseInverse());
    //timer[21] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    //st1 = clock();
    //temp = emPars.C[g].transpose() * temp;
    //timer[22] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    //st1 = clock();
    //temp /= temp.norm();
    //timer[23] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    st1 = clock();
    //Eigen::SparseMatrix<double> CC = emPars.C[g].transpose() * emPars.C[g];
    timer[19] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    st1 = clock();
    //Eigen::SparseMatrix<double> CM = emPars.C[g].transpose() * emPars.Mu[g];
    timer[20] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    st1 = clock();
    Eigen::SparseMatrix<double> CCX = CC * sX;
    CCX.prune(epsZero, 0.01);
    CCX.makeCompressed();
    timer[21] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    st1 = clock();
    Eigen::MatrixXd temp = CM.cwiseProduct(CCX.cwiseInverse());
    timer[22] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    st1 = clock();
    temp /= temp.norm();
    timer[23] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    QGrad_t += (double)(clock() - st) / CLOCKS_PER_SEC;
    QGrad_n ++;
    return temp;
}
