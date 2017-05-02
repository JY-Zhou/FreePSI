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
        string expName, int readLength, int THREAD) {
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
    time_t st, ed;
    time(&st);
    KmerHash kmerHasher(K, genomePath, exonBoundaryPath, readPath, readLength);
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
    cout << "FreePSI v0.1" << endl;
    cout << endl;
    cout << "Usage:" << endl;
    cout << "freePSI [options] -g <GENOME_DIR> -a <EXON_BND_ANNOT> [-r <KMER_COUNT> | -1 <KMER_COUNT> -2 <KMER_COUNT>] -o <OUTPUT>" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "  -g <GENOME_DIR>                  The directory containing the reference genome of each chromosome (.fasta format)" << endl;
    cout << "  -a <EXON_BND_ANNOT>              The annotation of exon boundary (.bed format)" << endl;
    cout << "  -r <KMER_COUNT>                  The k-mer count of single-end reads  produced by Jellyfish (.fasta format)" << endl;
    cout << "  -1 <KMER_COUNT> -2 <KMER_COUNT>  The k-mer count of paired-end reads produced by Jellyfish (.fasta format)" << endl;
    cout << "  -o <OUTPUT>                      The result of PSI values" << endl;
    cout << "  -k [Integer]                     The length of k-mer (Default: 27)" << endl;
    cout << "  -p [Thread number]               The thread numbers. (Default: 1)" << endl;
    cout << "  -h                               Help information" << endl;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);  
    std::cin.tie(0);

    string genomePath = "";
    string exonBoundaryPath = "";
    vector <string> readPath;
    int parameterCode = 0;
    int K = 27;
    int THREAD = 1;
    int readLength = 30;
    string expName = "";

    for(int i = 1; i < argc; i ++) {
        if(strcmp(argv[i], "-g") == 0) {
            genomePath = argv[++i];
            parameterCode |= 0x0001;
        } else if(strcmp(argv[i], "-a") == 0) {
            exonBoundaryPath = argv[++i];
            parameterCode |= 0x0002;
        } else if(strcmp(argv[i], "-r") == 0) {
            readPath.push_back(argv[++i]);
            parameterCode |= 0x0004;
        } else if(strcmp(argv[i], "-1") == 0) {
            readPath.push_back(argv[++i]);
            parameterCode |= 0x0008;
        } else if(strcmp(argv[i], "-2") == 0) {
            readPath.push_back(argv[++i]);
            parameterCode |= 0x0010;
        } else if(strcmp(argv[i], "-k") == 0) {
            K = boost::lexical_cast<int>(argv[++i]);
            parameterCode |= 0x0020;
        } else if(strcmp(argv[i], "-l") == 0) {
            readLength = boost::lexical_cast<int>(argv[++i]);
            parameterCode |= 0x0040;
        } else if(strcmp(argv[i], "-o") == 0) {
            expName = argv[++i];
            parameterCode |= 0x0080;
        } else if(strcmp(argv[i], "-p") == 0) {
            THREAD = boost::lexical_cast<int>(argv[++i]);
            parameterCode |= 0x0100;
        } else if(strcmp(argv[i], "-h") == 0) {
            showHelp();
            return 0;
        }
    }

    if(parameterCode == 0) {
        showHelp();
        exit(-1);
    } else if((parameterCode & 0x0001) == 0) {
        cerr << "\nError: Require reference genome (missing option -g)!" << endl;
        cerr << "\nPlease use -h option to show help information.\n" << endl;
        exit(-1);
    } else if((parameterCode & 0x0002) == 0) {
        cerr << "\nError: Require exon boundary annotation (missing option -b)!" << endl;
        cerr << "\nPlease use -h option to show help information.\n" << endl;
        exit(-1);
    } else if((parameterCode & 0x001c) == 0) {
        cerr << "\nError: Require k-mer count (missing option -r or -1 -2)!" << endl;
        cerr << "\nPlease use -h option to show help information.\n" << endl;
        exit(-1);
    } else if((parameterCode & 0x001c) != 0x0018 && (parameterCode & 0x001c) != 0x0004) {
        cerr << "\nError: K-mer count input is confusing (-r or -1 -2)" << endl;
        cerr << "\nPlease use -h option to show help information.\n" << endl;
        exit(-1);
    }

    long st = clock();
    time_t n_st, n_ed;
    time(&n_st);
    work(genomePath, exonBoundaryPath, readPath, K, expName, readLength, THREAD);
    time(&n_ed);


    cout << "CPU Time elapsed: " << (clock() - st) / CLOCKS_PER_SEC << "s. "<< endl;
    cout << "Natural Time elapsed: " << difftime(n_ed, n_st) << "s. " << endl;

    return 0;
}
