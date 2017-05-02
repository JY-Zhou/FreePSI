#include "BiasCorrection.h"

BiasCorrection::BiasCorrection() {}

BiasCorrection::BiasCorrection(KmerHash& kmerHasher, EMAlgorithm& EMSolution) {
    geneBoundary = kmerHasher.geneBoundary;
    genomePath = kmerHasher.genomePath;
    NG = EMSolution.NG;
    NE = EMSolution.NE;
    NX = EMSolution.NX;
    NC = 1 + 2 + 4 * 4;
    readLength = kmerHasher.readLength;

    for(int g = 0; g < NG; g ++) {
        X.push_back(Eigen::MatrixXd::Zero(NX[g], NC));
        P.push_back(Eigen::MatrixXd::Zero(NX[g], 0));
        Beta.push_back(Eigen::MatrixXd::Zero(NX[g], 1));
        Y.push_back(Eigen::MatrixXd::Zero(NX[g], 1));
    }
    initialX();
    initialY(EMSolution);
}

BiasCorrection::~BiasCorrection() {}

inline int BiasCorrection::parseBP(char c) {
    if(c == 'A' or c == 'a') return 0;
    if(c == 'C' or c == 'c') return 1;
    if(c == 'G' or c == 'g') return 2;
    if(c == 'T' or c == 't') return 3;
    return -1;
}

void BiasCorrection::initialX() {
    ifstream genomeFile(genomePath.c_str(), ios::in);
    string geneSeq = "";
    string line;
    while(getline(genomeFile, line)) {
        if(line[0] != '>') {
            geneSeq += boost::trim_copy(line);
        }
    }
    genomeFile.close();

    long long mask = (1L << 4) - 1;

    for(int g = 0; g < NG; g ++) {
        int row = 0;
        
        for(int e = 0; e < NE[g]; e ++) {
            row = e;

            int gcContent = 0;
            int totId = 0;
            X[g](row, 0) = 1.0;

            int exonst = geneBoundary[g][e].first;
            int exoned = geneBoundary[g][e].second + 1;

            X[g](row, 1) = exoned - exonst;

            int p = 0;
            long long id = 0;

            for(; p < 1; p ++) {
                int val = parseBP(geneSeq[exonst + p]);
                id = (id << 2) + val;
                if(val == 1 || val == 2) {
                    gcContent ++;
                }
            }

            for(; exonst + p < exoned; p ++) {
                int val = parseBP(geneSeq[exonst + p]);
                id = ((id << 2) + val) & mask;
                X[g](row, 3 + id) ++;
                totId ++;
                if(val == 1 || val == 2) {
                    gcContent ++;
                }
            }

            X[g](row, 2) = gcContent / X[g](row, 1);
            for(int c = 0; c < 16; c ++) {
                X[g](row, 3 + c) /= totId;
            }
        }

        for(int ei = 0; ei < NE[g]; ei ++) {
            for(int ej = ei + 1; ej < NE[g]; ej ++) {
                row = (2 * NE[g] - ei) * (ei + 1) / 2 + ej - ei - 1;

                int gcContent = 0;
                int totId = 0;

                X[g](row, 0) = 1.0;
                
                int exonedi = geneBoundary[g][ei].second + 1;
                int exonstj = geneBoundary[g][ej].first;

                X[g](row, 1) = 2 * readLength -  2;

                int p = 0;
                long long id = 0;

                for(; p < 1; p ++) {
                    int val = parseBP(geneSeq[exonedi - readLength + 1 + p]);
                    id = (id << 2) + val;
                    if(val == 1 || val == 2) {
                        gcContent ++;
                    }
                }

                for(; p < readLength - 1; p ++) {
                    int val = parseBP(geneSeq[exonedi - readLength + 1 + p]);
                    id = ((id << 2) + val) & mask;
                    X[g](row, 3 + id) ++;
                    totId ++;
                    if(val == 1 || val == 2) {
                        gcContent ++;
                    }
                }

                for(p = 0; p < readLength - 1; p ++) {
                    int val = parseBP(geneSeq[exonstj + p]);
                    id = ((id << 2) + val) & mask;
                    X[g](row, 3 + id) ++;
                    totId ++;
                    if(val == 1 || val == 2) {
                        gcContent ++;
                    }
                }
                 
                X[g](row, 2) = gcContent / X[g](row, 1);
                for(int c = 0; c < 16; c ++) {
                    X[g](row, 3 + c) /= totId;
                }
            }
        }

        for(int i = 0; i < X[g].rows(); i ++) {
            X[g](i, 1) /= X[g].col(1).sum();
        }
    }
    return;
}

void BiasCorrection::initialY(EMAlgorithm& EMSolution) {
    for(int g = 0; g < NG; g ++) {
        Y[g] = EMSolution.X[g];
    }
    return;
}

void BiasCorrection::PCA(int g) {
    Eigen::MatrixXd centered = X[g].rowwise() - X[g].colwise().mean();
    Eigen::MatrixXd cov = (centered.adjoint() * centered) / (double)(X[g].rows() - 1);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
    
    // cout << eig.eigenvalues() << endl;
    // cout << eig.eigenvalues().sum() << endl;
    // cout << eig.eigenvalues().rows() << "," << eig.eigenvalues().cols() << endl;
    // cout << eig.eigenvectors().rows() << "," << eig.eigenvectors().cols() << endl;
    
    Eigen::MatrixXd transfer = Eigen::MatrixXd(NC, 0);

    double totVar = 0.0;
    for(int i = eig.eigenvalues().rows() - 1; i >= 0; i --) {
        totVar += eig.eigenvalues()(i, 0);

        transfer.conservativeResize(transfer.rows(), transfer.cols()+1);
        transfer.col(transfer.cols()-1) = eig.eigenvectors().col(i);

        if(totVar / eig.eigenvalues().sum() > 0.95) {
            break;
        }
    }
    P[g] = X[g] * transfer;

    return;
}

void BiasCorrection::regression(int g) {
    cout << "PC" << g << " --> " << P[g].cols() << endl;

    Beta[g] = P[g].jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Y[g]);
    Y[g] = P[g] * Beta[g];

    double tot = Y[g].sum();
    for(int i = 0; i < Y[g].rows(); i ++) {
        Y[g](i, 0) /= tot;
    }
    return;
}

void BiasCorrection::work() {
    cout << "=== Correcting bias.. === " << endl;
    for(int g = 0; g < NG; g ++) {
        if(NE[g] > 1) {
            PCA(g);
            regression(g);
        }
    }
    return;
}
