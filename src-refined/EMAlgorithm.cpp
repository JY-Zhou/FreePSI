#include <cstdlib>
#include "EMAlgorithm.h"

EMAlgorithm::EMAlgorithm(KmerHash& kmerHasher, int nThread): kmerHasher(kmerHasher) {
    cout << "\n\n\n### Start to initialize parameters ..." << endl;

    K = kmerHasher.K;
    readLength = kmerHasher.readLength;
    NW = kmerHasher.NW;
    NG = kmerHasher.NG;
    NE.clear();
    this -> nThread = nThread;
    for(int g = 0; g < NG; g ++) {
        NE.push_back(kmerHasher.NE[g]);
    }

    procSegment = NG / 50;

    time_t st, ed;
    time(&st);
    initialIndices();
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    initialCoefficients(kmerHasher);
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    initialVariables();
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    initialConstraints();
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    correctInfeasibleVariables();
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    cout << "\n### Finish initializing parameters !" << endl;
}

EMAlgorithm::~EMAlgorithm() {}

inline tuple <int, int, int> EMAlgorithm::decryptExonId(long long exonId) {
    return make_tuple(exonId >> 20, exonId >> 10 & ((1<<10)-1), exonId & ((1<<10)-1));
}

void EMAlgorithm::initialIndices() {
    cout << "\n*** Start to initialize variable indices ..." << endl;
    int proc = 0;
    NX.clear();
    for(int g = 0; g < NG; g ++) {
        if(g > proc) {
            cout << ">" << flush;
            proc += procSegment;
        }
        NX.push_back(NE[g] * (NE[g] + 1) / 2);
    }
    cout << " OK!\n*** Finish initializing variable indices!" << endl;
}

void EMAlgorithm::initialCoefficients(KmerHash& kmerHasher) { 
    cout << "\n*** Start to initialize coefficients ..." << endl;

    cout << "+++ Start to initialize effective length ..." << endl;
    int proc = 0;
    L.clear();
    for(int g = 0; g < NG; g ++) {
        if(g > proc) {
            cout << ">" << flush;
            proc += procSegment;
        }
        L.push_back(Eigen::MatrixXd(1, NX[g]));
        int col = 0;
        for(int e = 0; e < NE[g]; e ++) {
            int order_e = -1;
            if(kmerHasher.geneStrand[g] == "-") {
                order_e = NE[g] - e - 1;
            } else if(kmerHasher.geneStrand[g] == "+") {
                order_e = e;
            } else {
                cout << "Error: Illegal strand info!" << endl;
                exit(-1);
            }

            int st = kmerHasher.geneBoundary[g][order_e].first;
            int ed = kmerHasher.geneBoundary[g][order_e].second;

            L[g](0, col) = ed - st - readLength + 1;
            if(L[g](0, col) <= 0) {
                cout << "Error: Short exon still exists!" << endl;
                exit(-1);
            }
            ++ col;
        }
        for(int ei = 0; ei < NE[g]; ei ++) {
            for(int ej = ei + 1; ej < NE[g]; ej ++) {
                L[g](0, col) = readLength - 1;
                ++ col;
            }
        }
    }
    cout << " OK!\n+++ Finish initializing effective length!" << endl;

    cout << "+++ Start to initialize contribution matrix and kmer count vector ..." << endl;
    C.clear();
    for(int g = 0; g < NG; g ++) {
        C.push_back(Eigen::SparseMatrix<double>(NW, NX[g]));
    }
    W.resize(NW, 1);
    int tempProcSegment = kmerHasher.kmerCount.size() / NSEGMENT;
    proc = 0;
    int row = 0;
    for(auto it: kmerHasher.kmerTable) {
        if(row > proc) {
            cout << ">" << flush;
            proc += tempProcSegment;
        }
        long long kmerId = it.first;
        W(row, 0) = kmerHasher.kmerCount[kmerId];
        for(auto exon: it.second) {
            int g, ei, ej;
            tie(g, ei, ej) = decryptExonId(exon.first);
            int col = 0;
            if(ei == ej) {
                col = ei;
            } else {
                col = (2 * NE[g] - ei) * (ei + 1) / 2 + ej - ei - 1;
            }
            C[g].insert(row, col) = exon.second;
        }
        ++ row;
    }
    for(int g = 0; g < NG; g ++) {
        C[g].prune(EPS, 0.01);
        C[g].makeCompressed();
    }
    cout << " OK!\n+++ Finish initializing contribution matrix and kmer count vector! " << endl;
    cout << "*** Finish initializing coefficients!" << endl;
}

void EMAlgorithm::initialVariables() {
    cout << "\n*** Start to set initial values for variables (Gamma, Theta and Mu)..." << endl;

    cout << "+++ Start to allocate kmer count as initial value ..." << endl;
    //alloc mem for variables
    Z.resize(NG, 1);
    X.clear();
    for(int g = 0; g < NG; g ++) {
        X.push_back(Eigen::SparseMatrix<double>(NX[g], 1));
    }
    Mu = vector<Eigen::SparseMatrix<double> >(NG, Eigen::SparseMatrix<double>(NW, 1));

    //update gene state
    int nManyExons = 0;
    isAvaliableZ = vector<int>(NG, 0);
    for(int g = 0; g < NG; g ++) {
        if(NE[g] >= exonNumBound) {
            nManyExons ++;
            isAvaliableZ[g] = 2;
        }
    }

    //init variables
    int nAverageShare = 0;
    int nUniqueManyExons = 0;
    int proc = 0, row = 0, tempProcSegment = kmerHasher.kmerCount.size() / NSEGMENT;
    for(auto it: kmerHasher.kmerTable) {
        if(row > proc) {
            cout << ">" << flush;
            proc += tempProcSegment;
        }

        nAverageShare += it.second.size();
        long long kmerId = it.first;
        double kmerNum = kmerHasher.kmerCount[kmerId];

        double sumContribution = 0.0;
        for(auto exon: it.second) {
            int g, ei, ej;
            tie(g, ei, ej) = decryptExonId(exon.first);
            int j = 0;
            if(ei == ej) {
                j = ei;
            } else {
                j = (2 * NE[g] - ei) * (ei + 1) / 2 + ej - ei - 1;
            }
            if(isAvaliableZ[g] != 2) { 
                sumContribution += C[g].coeffRef(row, j);
            }
        }

        for(auto exon: it.second) {
            int g, ei, ej;
            tie(g, ei, ej) = decryptExonId(exon.first);
            int j = 0;
            if(ei == ej) {
                j = ei;
            } else {
                j = (2 * NE[g] - ei) * (ei + 1) / 2 + ej - ei - 1;
            }

            if(isAvaliableZ[g] != 2) {
                isAvaliableZ[g] = 1;
                Z.coeffRef(g, 0) += (kmerNum / sumContribution) * C[g].coeffRef(row, j);
                X[g].coeffRef(j, 0) += (kmerNum / sumContribution) * C[g].coeffRef(row, j);
            } else {
                if(it.second.size() == 1) {
                    nUniqueManyExons ++;
                }
            }
        }
        row ++;
    }

    Z.prune(EPS, 0.01);
    Z.makeCompressed();
    for(int g = 0; g < NG; g ++) {
        X[g].prune(EPS, 0.01);
        X[g].makeCompressed();
    }

    cout << " OK!\n+++ Finish allocating kmer count as initial value ..." << endl;

    int nTotalX = 0;
    int nZeroExon = 0;
    int nTotalExon = 0;
    int nZeroX = 0;
    int nZeroZ = NG - Z.nonZeros();
    for(int g = 0; g < NG; g ++) {
        nTotalX += NX[g];
        nZeroX += NX[g] - X[g].nonZeros();
        if(isAvaliableZ[g] != 2) {
            for(int e = 0; e < NE[g]; e ++) {
                nTotalExon ++;
                if(X[g].coeffRef(e, 0) < EPS) {
                    nZeroExon ++;
                }
            }
        }
    }
    cout << "=== On average, " << (double)nAverageShare / NW << " exons or juctions share one kmer." << endl;
    cout << "=== " << (double)nUniqueManyExons / NW * 100 << "\% kmers are unique to exons in large exon number gene." << endl;
    cout << "=== " << (double)nZeroZ / NG * 100 << "\% elements of Z (" << NG << " in total) are zero." << endl; 
    cout << "=== " << (double)nZeroX / nTotalX * 100 << "\% elements of X (" << nTotalX << " in total) are zero." << endl;
    cout << "=== " << (double)nZeroExon << " elements of X (" << nTotalExon << " in total) are zero." << endl;
    cout << "=== " << (double)nManyExons / NG * 100 << "\% genes (" << nManyExons << " in total) contain more than " << exonNumBound << " exons"  << endl;

    cout << "+++ Start to do normalization ..." << endl;
    proc = 0;
    Z /= Z.sum();
    for(int g = 0; g < NG; g ++) {
        if(g > proc) {
            cout << ">" << flush;
            proc += procSegment;
        }
        if(isAvaliableZ[g] == 1) {
            X[g] /= X[g].sum();
        }
    }
    cout << " OK!\n+++ Finish doing normalization!" << endl;
    cout << "*** Finish setting initial values for variables (Gamma, Alpha and Mu)!" << endl;
}

void EMAlgorithm::initialConstraints() {
    cout << "\n*** Start to initialize constraint matrix ..." << endl;
    int proc = 0;

    NA.clear();
    for(int g = 0; g < NG; g ++) {
        if(NE[g] > 1) {
            NA.push_back(NE[g] + NX[g] - 2);
        } else {
            NA.push_back(NE[g]);
        }
    }

    A.clear();

    for(int g = 0; g < NG; g ++) {
        if(g > proc) {
            cout << ">" << flush;
            proc += procSegment;
        }

        A.push_back(Eigen::SparseMatrix<double>(NA[g], NX[g]));

        if(isAvaliableZ[g] == 1) {
            if(NA[g] == 1) {
                A[g].insert(0, 0) = 1 / L[g](0, 0);
                continue;
            }

            int row = 0;

            int NRightJunc = NE[g] - 1;
            int l = NE[g];
            int r = l + NE[g] - 1;
            for(; row < NRightJunc ;row ++) {
                A[g].insert(row, row) = 1 / L[g](0, row);
                for(int i = l; i < r; i ++) {
                    A[g].insert(row, i) = -1 / L[g](0, i);
                }
                l = r;
                r = r + NE[g] - (row + 2);
            }

            int NLeftJunc = NRightJunc + NE[g] - 1;
            for(; row < NLeftJunc ; row ++) {
                A[g].insert(row, row - NRightJunc + 1) = 1 / L[g](0, row - NRightJunc + 1);
                l = NE[g];
                r = row - NRightJunc;
                for(int i = 1; r >= 0; i ++, r --) {
                    A[g].insert(row, l + r) = -1 / L[g](0, l + r);
                    l += NE[g] - i;
                }
            }

            for(int i = 0; i < NX[g] - NE[g]; i ++, row ++) {
                A[g].insert(row, NE[g] + i) = 1;
            }
            A[g].prune(EPS, 0.01);
            A[g].makeCompressed();
        }
    }
    cout << "OK! \n*** Finish initializing constraint matrix!" << endl;
}

void EMAlgorithm::correctInfeasibleX(int g) {
    X[g] /= X[g].sum();

    int t = 0;
    while(true) {
        X[g].prune(EPS, 0.01);
        bool ok = true;
        Eigen::SparseMatrix<double> checkIEQ = A[g] * X[g];
        for(int i = 0; i < checkIEQ.rows(); i ++) {
            if(checkIEQ.coeffRef(i, 0) < - EPS) {
                ok = false;
                int xid = -1;
                if(i < NE[g] - 1) xid = i;
                else if(i < 2*NE[g] - 2) xid = i - NE[g] + 2; 

                if(fabs(X[g].coeffRef(xid, 0)) < EPS) {
                    for(int j = 0; j < NX[g]; j ++) {
                        if(fabs(A[g].coeffRef(i, j)) > EPS) {
                            X[g].coeffRef(j, 0) = 0;
                        }
                    }
                }
            }
        }
        if(ok) break;
        X[g].prune(EPS, 0.01);
        for(int i = 0; i < NX[g]; i ++) {
            if(i < NE[g]) {
                X[g].coeffRef(i, 0) *= 10;
            } 
        }
    }

    X[g] /= X[g].sum();
    X[g].prune(EPS, 0.01);
    X[g].makeCompressed();
}

void EMAlgorithm::correctInfeasibleVariables() {
    cout << "\n*** Start to make Alpha feasible..." << endl;
    int proc = 0;
    int nShrinkage = 0;
    int nOriginal = 0;
    int nEmpty = 0;
    for(int g = 0; g < NG; g ++) {
        if(g > proc) {
            cout << ">" << flush;
            proc += procSegment;
        }

        if(isAvaliableZ[g] == 1) {
            int nzz = X[g].nonZeros();
            correctInfeasibleX(g);
            int projectedNZZ = X[g].nonZeros();

            if(projectedNZZ > nzz) {
                cout << "Error: variable shrinkage wont increase nonzeros!" << endl;
                exit(-1);
            }

            nShrinkage += nzz - projectedNZZ;
            nOriginal += nzz;

            if(projectedNZZ == 0) {
                nEmpty ++;
                isAvaliableZ[g] = 0;
                continue;
            }

            //check
            Eigen::SparseMatrix<double> ieq = A[g] * X[g];
            for(int i = 0; i < ieq.rows(); i ++) {
                if(ieq.coeffRef(i, 0) < - EPS ) {
                    cout << "Error: Still infeasible point! (Violate inequalities)" << endl;
                    cout << ieq.coeffRef(i, 0) << " in gene " << g << " constraint " << i << endl;
                }
            }
            Eigen::MatrixXd eq = Eigen::MatrixXd::Ones(1, NX[g]) * X[g];
            if(fabs(eq(0, 0) - 1) > EPS ) {
                cout << "Error: Still infeasible point! (Violate equalities)" << endl;
                cout << eq(0, 0) << " in gene " << g << endl;
            }

        }
    }


    cout << "OK! \n+++ Finish making Alpha feasible!" << endl;
    cout << "=== Prune " << (double)nShrinkage / nOriginal * 100 << "\% infeaible variables (infeasible junctions) ..." << endl; 
    cout << "=== " << nEmpty << " new empty genes are found" << endl; 
    cout << "=== Now there are " << nOriginal - nShrinkage << " variables. (Before " << nOriginal << ")" << endl;
            
    cout << "*** Finish initializing CGPA optimizer and making Alpha feasible..." << endl;
}

#pragma omp declare reduction (\
        +: \
        Eigen::VectorXd: \
        omp_out = omp_out + omp_in) \
initializer(omp_priv = Eigen::VectorXd::Zero(omp_orig.size()))

//bool DEBUG = false;

void EMAlgorithm::eStep() {
    cout << "+++ E-step ..." << endl;
    int proc = 1;
    vector <double> timer = vector <double> (10, 0);

    clock_t st = clock();
    Eigen::VectorXd sumMu = Eigen::VectorXd::Zero(NW);

    clock_t st1 = clock();
//#ifdef EIGEN_DONT_PARALLELIZE
//#pragma omp parallel for num_threads(1)
//#endif
    for(int g = 0; g < NG; g ++) {
//#ifdef EIGEN_DONT_PARALLELIZE
//#pragma omp critical
//#endif
        {
            if(g > proc) {
                cout << ">" << flush;
                proc += procSegment * 2;
            }
        }

        if(isAvaliableZ[g] == 1) {
            Mu[g].setZero();
            Mu[g].makeCompressed();
            for(Eigen::SparseMatrix<double>::InnerIterator it(X[g], 0); it; ++ it) {
                Mu[g] +=  C[g].col(it.row()) * it.value() * Z.coeffRef(g, 0);
            }

            Mu[g].prune(EPS, 0.01);
            Mu[g].makeCompressed();
        }
    }

#pragma omp parallel for num_threads(nThread) reduction(+:sumMu)
    for(int g = 0; g < NG; g ++) {
        if(isAvaliableZ[g] == 1) {
            sumMu += Mu[g];
        }
    }
    timer[1] += (double)(clock() - st1) / CLOCKS_PER_SEC;

    //if(DEBUG) exit(-1);
        
    clock_t st2 = clock();
    timer[2] += (double)(clock() - st2) / CLOCKS_PER_SEC;

    proc = 0;
    for(int g = 0; g < NG; g ++) {
        if(g > proc) {
            cout << ">" << flush;
            proc += procSegment * 2;
        }

        if(isAvaliableZ[g] == 1) {
            clock_t st3 = clock();
            Mu[g] = Mu[g].cwiseProduct(sumMu.cwiseInverse());
            timer[3] += (double)(clock() - st3) / CLOCKS_PER_SEC;

            clock_t st4 = clock();
            Mu[g] = Mu[g].cwiseProduct(W);
            timer[4] += (double)(clock() - st4) / CLOCKS_PER_SEC;

            clock_t st5 = clock();
            Mu[g].prune(EPS, 0.01);
            Mu[g].makeCompressed();
            timer[5] += (double)(clock() - st5) / CLOCKS_PER_SEC;
        }
    }
    timer[0] += (double)(clock() - st) / CLOCKS_PER_SEC;
    cout << " OK!" << endl;
}

void EMAlgorithm::mStep(double eps, int iter) {
    cout << "Thread number = " << nThread << endl;
    cout << "+++ M-step ..." << endl;
    for(int g = 0; g < NG; g ++) {
        if(isAvaliableZ[g] == 1) {
            Z.coeffRef(g, 0) = Mu[g].sum();
            if(Z.coeffRef(g, 0) < EPS) {
                isAvaliableZ[g] = 0;
            }
        }
    }
    Z /= Z.sum();

    for(int i = 0; i < ConjugateGradientProjection::timer.size(); i ++) {
        ConjugateGradientProjection::timer[i] = 0;
    }
    cout << "Eigen threads = " << Eigen::nbThreads() << endl;

    int proc = 0;
#ifdef EIGEN_DONT_PARALLELIZE
#pragma omp parallel for num_threads(nThread)
#endif
    for(int g = 0; g < NG; g ++) {
#ifdef EIGEN_DONT_PARALLELIZE
#pragma omp critical
#endif
        {
            if(g == 0) cout << "Thread num = " << omp_get_num_threads() << endl;
            if(g > proc) {
                cout << ">" << flush;
                proc += procSegment;
            }
        }

        if(isAvaliableZ[g] == 1) {
            int bakNA = NA[g];
            int bakNX = NX[g];
            Eigen::SparseMatrix<double> bakA = A[g];
            Eigen::SparseMatrix<double> bakC = C[g];
            Eigen::SparseMatrix<double> bakX = X[g];
            vector<int> nonzeroXOriginalIndex;

            nonzeroXOriginalIndex.clear();
            bakX.prune(EPS, 0.01), bakX.makeCompressed();
            for(int i = 0; i < bakNX; i ++) {
                if(bakX.coeffRef(i, 0) > EPS * 1e-2) {
                    nonzeroXOriginalIndex.push_back(i);
                }
            }

            NX[g] = nonzeroXOriginalIndex.size();

            A[g].resize(NA[g], NX[g]);
            C[g].resize(NW, NX[g]);
            X[g].resize(NX[g], 1);
            for(int i = 0; i < nonzeroXOriginalIndex.size(); i ++) {
                int j = nonzeroXOriginalIndex[i];
                A[g].col(i) = bakA.col(j);
                C[g].col(i) = bakC.col(j);
                X[g].insert(i, 0) = bakX.coeffRef(j, 0);
            }
            A[g].prune(EPS, 0.01), A[g].makeCompressed();
            C[g].prune(EPS, 0.01), C[g].makeCompressed();
            X[g].prune(EPS, 0.01), X[g].makeCompressed();

            Eigen::SparseMatrix<double, Eigen::RowMajor> tempA = A[g];
            Eigen::SparseMatrix<double, Eigen::RowMajor> newA(0, NX[g]);
            int nANZRows = 0;
            for(int i = 0; i < bakNA; i ++) {
                if(tempA.row(i).norm() > EPS * 1e-2) {
                    newA.conservativeResize(newA.rows() + 1, NX[g]);
                    newA.bottomRows<1>() = tempA.row(i);
                    nANZRows ++;
                }
            }
            NA[g] = nANZRows;
            A[g] = Eigen::SparseMatrix<double, Eigen::ColMajor>(newA);


            //Go to solver!
            Eigen::SparseMatrix<double> newX = X[g];
            if(X[g].nonZeros() > 0) {
                ConjugateGradientProjection optimizer(*this, g);
                newX = Eigen::MatrixXd(optimizer.optimizeQ(eps, iter)).sparseView();
            } else {
                cout << "Error!: but go on, after compact ,X is empty!" << endl;
                cout << "Gene = " << g << "NX = " << NX[g] << " bakNX = " << bakNX << endl;
                exit(-1);
            }

            newX.prune(EPS, 0.01), newX.makeCompressed();
            NX[g] = bakNX;
            NA[g] = bakNA;
            A[g] = bakA;
            C[g] = bakC;

            int nNonZeroX = 0;
            X[g].resize(bakNX, 1);
            for(int i = 0; i < nonzeroXOriginalIndex.size(); i ++) {
                int j = nonzeroXOriginalIndex[i];
                if(newX.coeffRef(i, 0) > EPS * 1e-2) {
                    nNonZeroX ++;
                    X[g].insert(j, 0) = newX.coeffRef(i, 0);
                }
            }
            A[g].prune(EPS, 0.01), A[g].makeCompressed();
            C[g].prune(EPS, 0.01), C[g].makeCompressed();
            X[g].prune(EPS, 0.01), X[g].makeCompressed();

            if(nNonZeroX == 0) {
                isAvaliableZ[g] = 0;
            }
        }
    }
    cout << " OK!" << endl;
}

double EMAlgorithm::computeLogLikelihood() {
    vector <double> timer = vector <double> (10, 0);
    clock_t st = clock();

    Eigen::VectorXd ZCX = Eigen::VectorXd::Zero(NW);


    timer[1] += (double)(clock() - st) / CLOCKS_PER_SEC;


    clock_t st2 = clock();

#ifdef EIGEN_DONT_PARALLELIZE
#pragma omp parallel for num_threads(nThread) reduction(+:ZCX)
#endif
    for(int g = 0; g < NG; g ++) {
        if(isAvaliableZ[g] == 1) {
            ZCX += C[g] * X[g] * Z.coeffRef(g, 0);
        }
    }
    
    timer[2] += (double)(clock() - st2) / CLOCKS_PER_SEC;


    clock_t st3 = clock();


    double ans = 0.0;
    for(int w = 0; w < NW; w ++) {
        if(ZCX(w, 0) > EPS) {
            ans += W(w, 0) * log(ZCX(w));
        }
    }
    
    timer[3] += (double)(clock() - st3) / CLOCKS_PER_SEC;

    timer[0] += (double)(clock() - st) / CLOCKS_PER_SEC;
    cout << " OK!" << endl;

    return ans;
}

void EMAlgorithm::work(string chrName) {
    cout << "\n\n\n### Start to estimate parameters ... " << endl;
    cout << "\n*** Start to run EM algorithm ..." << endl;
    cout << "*** M-step mode:" << endl;
#ifdef EIGEN_DONT_PARALLELIZE
    cout << "--- Run M-step concurrently" << endl;
#else
    cout << "--- Dont run M-step concurrently" << endl;
#endif

    double cpu_EStep = 0.0, cpu_MStep = 0.0, cpu_Likelihood = 0.0;
    double nat_EStep = 0.0, nat_MStep = 0.0, nat_Likelihood = 0.0;

    int maxIter = 50;

    time_t st, ed;
    clock_t st_Likelihood = clock();
    time(&st);

    preLikelihood = computeLogLikelihood();

    cpu_Likelihood += (double)(clock() - st_Likelihood) / CLOCKS_PER_SEC;
    time(&ed);
    nat_Likelihood += difftime(ed, st);

    Eigen::SparseMatrix <double> preZ = Z;
    cout << "=== Initial log-likelihood is " << preLikelihood << endl;
    for(int it = 1; it <= maxIter; it ++) {
        cout << "\n+++ " << it << " iteration processed..." << endl;
        
        clock_t st_EStep = clock();
        time(&st);

        //if(it == 4) DEBUG = true;

        eStep();

        cpu_EStep += (double)(clock() - st_EStep) / CLOCKS_PER_SEC;
        time(&ed);
        nat_EStep += difftime(ed, st);
        cout << "Elasped time " << difftime(ed, st) << "s. " << endl;


        clock_t st_MStep = clock();
        time(&st);

        mStep(1e-12, 100);

        cpu_MStep += (double)(clock() - st_MStep) / CLOCKS_PER_SEC;
        time(&ed);
        nat_MStep += difftime(ed, st);
        cout << "Elasped time " << difftime(ed, st) << "s. " << endl;


        st_Likelihood = clock();
        time(&st);

        newLikelihood = computeLogLikelihood();

        cpu_Likelihood += (double)(clock() - st_Likelihood) / CLOCKS_PER_SEC;
        time(&ed);
        nat_Likelihood += difftime(ed, st);

        cout << "=== New log-likelihood is " << newLikelihood <<
            " (delta = " << newLikelihood - preLikelihood << ", " <<
            (newLikelihood - preLikelihood) / fabs(preLikelihood) * 100 << "% )." << endl;
        cout << "=== The L2-norm of difference of hidden variable is " << (preZ - Z).norm() << endl;
        if((newLikelihood - preLikelihood) / fabs(preLikelihood) < 1e-6) {
            for(int g = 0; g < NG; g ++) {
                if(isAvaliableZ[g] == 1 && Z.coeffRef(g, 0) < EPS) {
                    isAvaliableZ[g] = 0;
                }
            }
            break;
        } else {
            preLikelihood = newLikelihood;
            preZ = Z;
        }
    }
    cout << "*** Finish EM algorithm !" << endl;
    cout << "### Finish estimating parameters !" << endl;
    cout << "### Compute Log-Likelihood " << cpu_Likelihood << " (CPU), " <<
        nat_Likelihood << " (Natural)." << endl;
    cout << "### E-step " << cpu_EStep << " (CPU), " << nat_EStep << " (Natural)." << endl;
    cout << "### M-step " << cpu_MStep << " (CPU), " << nat_MStep << " (Natural)." << endl;
}

void EMAlgorithm::computePSI(string expName) {
    cout << "\n\n### Start to refine solution and compute PSI ... " << endl;
    PSI.clear();

    int proc = 0;
    for(int g = 0; g < NG; g ++) {
        if(g > proc) {
            cout << ">" << flush;
            proc += procSegment;
        }

        if(isAvaliableZ[g] == 1) {
            Eigen::SparseMatrix<double> tempPSI = X[g].cwiseProduct(L[g].transpose().cwiseInverse());
            double sumExon = 0.0, sumJunction = 0.0;
            for(int i = 0; i < NX[g]; i ++) {
                if(i < NE[g]) {
                    sumExon += tempPSI.coeffRef(i, 0);
                } else {
                    sumJunction += tempPSI.coeffRef(i, 0);
                }
            }
            tempPSI /= (sumExon - sumJunction);
            vector <double> temp;
            for(int e = 0; e < NE[g]; e ++) {
                temp.push_back(tempPSI.coeffRef(e, 0));
            }
            if(kmerHasher.geneStrand[g] == "-") {
                reverse(temp.begin(), temp.end());
            } 
            PSI.push_back(temp);
        } else if(isAvaliableZ[g] == 2){
            vector <double> temp(NE[g], 0.0);
            if(kmerHasher.geneStrand[g] == "-") {
                reverse(temp.begin(), temp.end());
            } 
            PSI.push_back(temp);
        } else {
            vector <double> temp(NE[g], 0.0);
            if(kmerHasher.geneStrand[g] == "-") {
                reverse(temp.begin(), temp.end());
            } 
            PSI.push_back(temp);
        }
    }

    string psiFile = expName + "/psi_freePSI_raw.json";
    ofstream psiOutput(psiFile, ios::out);
    psiOutput << "[" << endl;
    for(int g = 0; g < NG; g ++) {
        psiOutput << " [" << endl;
        for(int e = 0; e < NE[g]; e ++) {
            if(e == NE[g] - 1) {
                psiOutput << "  " << PSI[g][e] << endl << " ]";
            } else {
                psiOutput << "  " << PSI[g][e] << "," << endl;
            }
        }
        if(g == NG - 1) {
            psiOutput << endl << "]" << endl;
        } else {
            psiOutput << "," << endl;
        }
    }
    psiOutput.flush();
    psiOutput.close();


    cout << " OK!\n### Finish refining solution and computing PSI ... " << endl;
}
