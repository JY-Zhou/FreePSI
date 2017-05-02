#include <iostream>
#include <algorithm>
#include <ctime>
#include <vector>
#include <iomanip>
#include <cstring>

#include <Eigen/Core>

#include <boost/lexical_cast.hpp>

#include "KmerHash.h"
#include "EMAlgorithm.h"
#include "BiasCorrection.h"

using namespace std;

KmerHash load(string dumpPath, int K, int readLength) {
    KmerHash kmerHasher;
    kmerHasher.K = K;
    kmerHasher.readLength = readLength;

    ifstream input(dumpPath.c_str(), ios::in);
    input >> kmerHasher.NW;
    input >> kmerHasher.NG;
    for(int g = 0; g < kmerHasher.NG; g ++) {
        int a;
        input >> a;
        kmerHasher.NE.push_back(a);
    }
    for(int w = 0; w < kmerHasher.NW; w ++) {
        long long a, b;
        input >> a >> b;
        kmerHasher.kmerCount[a] = b;
    }
    for(int w = 0; w < kmerHasher.NW; w ++) {
        long long a;
        int s;
        input >> a >> s;
        unordered_map<long long, double> temp;
        for(int i = 0; i < s; i ++) {
            long long b;
            double c;
            input >> b >> c;
            temp[b] = c;
        }
        kmerHasher.kmerTable[a] = temp;
    }
    for(int g = 0; g < kmerHasher.NG; g ++) {
        vector<pair<int, int> > temp;
        for(int e = 0; e < kmerHasher.NE[g]; e ++) {
            int a, b;
            input >> a >> b;
            temp.push_back(make_pair(a, b));
        }
        kmerHasher.geneBoundary.push_back(temp);
    }
    for(int g = 0; g < kmerHasher.NG; g ++) {
        string s;
        input >> s;
        kmerHasher.geneStrand.push_back(s);
    }
    kmerHasher.kmerCount.rehash(kmerHasher.kmerCount.size());
    kmerHasher.kmerTable.rehash(kmerHasher.kmerTable.size());
    return kmerHasher;
}

void save(string dumpPath, KmerHash& kmerHasher) {
    ofstream output(dumpPath.c_str(), ios::out);
    output << kmerHasher.NW << endl;
    output << kmerHasher.NG << endl;
    for(int i = 0; i < kmerHasher.NE.size(); i ++) {
        output << kmerHasher.NE[i] << " ";
    }
    output << endl;
    for(auto it: kmerHasher.kmerCount) {
        output << it.first << " " << it.second << endl;
    }
    for(auto it: kmerHasher.kmerTable) {
        output << it.first << " " << it.second.size() << endl;
        for(auto exonit: it.second) {
            output << exonit.first << " " << fixed << setprecision(9) <<exonit.second << " ";
        }
        output << endl;
    }
    for(int g = 0; g < kmerHasher.NG; g ++) {
        for(int e = 0; e < kmerHasher.NE[g]; e ++) {
            output << kmerHasher.geneBoundary[g][e].first << " " << kmerHasher.geneBoundary[g][e].second << " ";
        }
        output << endl;
    }
    for(int g = 0; g < kmerHasher.NG; g ++) {
        output << kmerHasher.geneStrand[g] << " ";
    }
    output << endl;
    output.flush();
    output.close();
}

void work(string genomePath, string exonBoundaryPath, vector <string> readPath, int K, 
        string expName, double lambda, int readLength, bool ifBiasCorrection, string predictExonPath, int THREAD) {
    time_t st, ed;
    time(&st);
    KmerHash kmerHasher(K, genomePath, exonBoundaryPath, readPath, readLength, predictExonPath);
    //save("./KmerHasherDump.out", kmerHasher);
    //KmerHash kmerHasher = load("./KmerHasherDump.out", K, readLength);
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    EMAlgorithm solver(kmerHasher, THREAD); 
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    solver.work(expName);
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;
    
    solver.computePSI(expName);
    cout << "### All Finished!" << endl;
}

void showHelp() {
    cout << "Allowed options" << endl;
    cout << endl;
    cout << "\t-g <Directory>\tThe directory contrains reference genome of each chromosome." << endl;
    cout << endl;
    cout << "\t-a <BED file>\tThe annotation of exons." << endl;
    cout << endl;
    cout << "\t-r <FASTQ file> | -1 <FASTQ file> -2 <FASTQ file>\tThe read sequences." << endl;
    cout << endl;
    cout << "\t-k <Integer>\tThe length of kmer when computing. (Default: 25)" << endl;
    cout << endl;
    cout << "\t-l <Integer>\tThe length of theoretical read when computing. (Default: 30)" << endl;
    cout << endl;
    cout << "\t-o <JSON file>\tThe result of the exon inclusion ratio." << endl;
    cout << endl;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);  
    std::cin.tie(0);

    cout << "\n*** Eigen mode: " << endl;
#ifdef EIGEN_DONT_PARALLELIZE 
    cout << "--- Not use intrinsic parallel in Eigen" << endl;
#else
    Eigen::initParallel();
    cout << "--- Use intrinsic parallel in Eigen" << endl;
#endif

#ifdef EIGEN_USE_MKL_ALL
    cout << "--- Use Intel MKL in Eigen" << endl;
#else
    cout << "--- Not use Intel MKL in Eigen" << endl;
#endif

    string genomePath = "";
    string exonBoundaryPath = "";
    vector <string> readPath;
    int K = 25;
    int THREAD = 16;
    int readLength = 30;
    string expName = "";
    double lambda;
    bool ifBiasCorrection = false;
    string predictExonPath = "";

    for(int i = 1; i < argc; i ++) {
        if(strcmp(argv[i], "-g") == 0) {
            genomePath = argv[++i];
        } else if(strcmp(argv[i], "-a") == 0) {
            exonBoundaryPath = argv[++i];
        } else if(strcmp(argv[i], "-r") == 0) {
            readPath.push_back(argv[++i]);
        } else if(strcmp(argv[i], "-1") == 0) {
            readPath.push_back(argv[++i]);
        } else if(strcmp(argv[i], "-2") == 0) {
            readPath.push_back(argv[++i]);
        } else if(strcmp(argv[i], "-k") == 0) {
            K = boost::lexical_cast<int>(argv[++i]);
        } else if(strcmp(argv[i], "-o") == 0) {
            expName = argv[++i];
        } else if(strcmp(argv[i], "-p") == 0) {
            THREAD = boost::lexical_cast<int>(argv[++i]);
        } else if(strcmp(argv[i], "-l") == 0) {
            readLength = boost::lexical_cast<int>(argv[++i]);
        } else if(strcmp(argv[i], "--pe") == 0) {
            predictExonPath = argv[++i];
        } else if(strcmp(argv[i], "--enable-bias-correction") == 0) {
            ifBiasCorrection = true;
        }
    }

    long st = clock();
    time_t n_st, n_ed;
    time(&n_st);
    work(genomePath, exonBoundaryPath, readPath, K, expName, lambda, readLength, ifBiasCorrection, predictExonPath, THREAD);
    time(&n_ed);

    cout << "CPU Time elapsed: " << (clock() - st) / CLOCKS_PER_SEC << "s. "<< endl;
    cout << "Natural Time elapsed: " << difftime(n_ed, n_st) << "s. " << endl;

    return 0;
}
