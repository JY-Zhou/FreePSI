#include "KmerHash.h"

KmerHash::KmerHash() {}

KmerHash::KmerHash(int K, string genomePath, string exonBoundaryPath, vector <string> allReadPaths,
        int readLength, string predictExonPath) {
    cout << "\n\n\n### Start to build theoretical and real kmer profile ..." << endl;

    srand(time(0));
    this -> genomePath = genomePath;
    this -> exonBoundaryPath = exonBoundaryPath;
    this -> allReadPaths = allReadPaths;

    strandSpec = true;
    BPQualBound = 20.0;
    AVGQualBound = 30.0;
    lowExpKmerBound = 0;

    this -> K = K;
    this -> readLength = readLength;
    mask = 0;
    for(int i = 0; i < K; i++) {
        mask = (mask << 2) | 3;
    }
    shiftLeftMost = 2 * K - 2; 

    NW = 0;
    NG = 0;
    NE.clear();
    kmerCount.clear();
    kmerTable.clear();
    geneBoundary.clear();
    geneStrand.clear();
    chromST.clear();
    chromName.clear();

    time_t st, ed;
    time(&st);
    readGenome(exonBoundaryPath, predictExonPath);
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;
    
    time(&st);
    buildKmerTable(genomePath);
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    readReads(allReadPaths);
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;

    time(&st);
    mergeKmerTable();
    time(&ed);
    cout << "Elasped time " << difftime(ed, st) << "s. " << endl;
    cout << "### Finish building theoretical and real kmer profile!" << endl;
}

KmerHash::~KmerHash() {}

inline long long KmerHash::parseBP(char c) {
    if(c == 'A' or c == 'a') return 0ll;
    if(c == 'C' or c == 'c') return 1ll;
    if(c == 'G' or c == 'g') return 2ll;
    if(c == 'T' or c == 't') return 3ll;
    return (long long)(rand() % 4);
}

inline long long KmerHash::parseAntiBP(char c) {
    return 3ll - parseBP(c);
}

inline long long KmerHash::encryptExonId(int g, int ei, int ej) {
    return (long long)g << 20 | (long long)ei << 10 | (long long) ej;
}

void KmerHash::readReads(vector <string> allReadPaths) {
    string ext = "";
    for(auto readPath : allReadPaths) {
        vector <string> fileInfo;
        boost::split(fileInfo, readPath, boost::is_any_of("."));
        string curExt = fileInfo.back();
        boost::algorithm::to_lower(curExt);
        if(ext.length() > 0 && ext != curExt) {
            cout << "Error: The format of reads files are inconsistent!" << endl;
            exit(-1);
        } else {
            ext = curExt;
        }
    }
    if(ext == "fq" || ext == "fastq") {
        readReadsFromFastq(allReadPaths);
    } else {
        readReadsFromJellyfish(allReadPaths);
        //cout << "Error: Currently only support .FASTQ format RNA-seq reads file!" << endl;
        //exit(-1);
    }
    return;
}

inline void KmerHash::updateKmerCount(long long id, long long number, long long & nFailKmer) {
    if(kmerTable.count(id) > 0) {
        if(kmerCount.count(id) > 0) {
            kmerCount[id] += number;
        } else {
            kmerCount[id] = number;
        }
    } else {
        nFailKmer += number;
    }
}

void KmerHash::readReadsFromJellyfish(vector <string> allReadPaths) {
    cout << "\n*** Start to load reads from Jellyfish..." << endl;
    long long nKmers = 0;
    long long nFilterOut = 0;
    long long nFailKmer = 0;
    long long proc = 0;

    long long number = 0;
    string kmer = "";

    kmerCount.reserve(kmerTable.size());

    for(int pairId = 0; pairId < allReadPaths.size(); pairId ++) {
        ifstream readsFile;
        char buf[1<<16];
        readsFile.rdbuf() -> pubsetbuf(buf, sizeof(buf));
        readsFile.open(allReadPaths[pairId].c_str(), ios::in);

        for(string line; getline(readsFile, line); ) {
            if(line[0] == '>') {
                number = boost::lexical_cast <long long> (boost::trim_copy(line).substr(1));
                getline(readsFile, line);
                kmer = boost::trim_copy(line);
            }

            if(kmer.find('N') == string::npos) {
                if(pairId == 0) {
                    long long id = 0;
                    for(int p = 0; p < K; p ++) {
                        char c = kmer[p];
                        id = ((id << 2) | parseBP(c)) & mask;
                    }
                    updateKmerCount(id, number, nFailKmer);
                } else if (pairId == 1){
                    long long antiId = 0;
                    for(int p = 0; p < K; p ++) {
                        char c = kmer[p];
                        antiId = ((antiId >> 2) | (parseAntiBP(c) << shiftLeftMost)) & mask;
                    }
                    updateKmerCount(antiId, number, nFailKmer);
                } else {
                    cout << "Error: The number of reads file will never larger than 2!" << endl;
                    exit(-1);
                }
            } else {
                nFilterOut += number;
            }

            nKmers += number;
            proc ++;
            if(proc % 1000000 == 0) {
                cout << "+++ " << proc / 1000000 << ",000,000 kmers loaded ..." << endl;
            }
        }

        readsFile.close();
    }

    cout << "=== Totally " << nKmers << " kmers loaded." << endl;
    cout << "=== Reserve " << nKmers - nFilterOut << " high quality Kmers. (Reserve " 
        << 100 - (double)nFilterOut / nKmers * 100 << "%)" << endl;
    cout << "=== " << nFailKmer << " kmers failed to match with theoretical kmer profile. (Reserve " 
        << 100 - (double)nFailKmer / nKmers * 100 << "%)" << endl;

    cout << "=== Finally collect " << kmerTable.size() << " different possible kmers on Genome."<< endl;
    cout << "=== Finally collect " << kmerCount.size() << " different kmers occur in read."<< endl;
    cout << "=== About " << (nKmers - - nFilterOut - nFailKmer) / kmerCount.size() << " occurrences of kmer on average." << endl;
    cout << "*** Finish loading kmers from reads!" << endl;
}

void KmerHash::readReadsFromFastq(vector <string> allReadPaths) {
    cout << "\n*** Start to load reads from raw FastQ file..." << endl;
    long long nRawReads = 0;
    long long nRealReads = 0;
    long long nFilterOut = 0;
    long long nFailKmer = 0;

    string identifier = "";
    string reads = "";
    string option = "";
    string score = "";

    kmerCount.reserve(kmerTable.size());

    //map <string, int> failset;
    for(int pairId = 0; pairId < allReadPaths.size(); pairId ++) {
        ifstream readsFile;
        char buf[1<<16];
        readsFile.rdbuf() -> pubsetbuf(buf, sizeof(buf));
        readsFile.open(allReadPaths[pairId].c_str(), ios::in);

        for(string line; getline(readsFile, line); ) {
            if(line[0] == '@') {
                identifier = boost::trim_copy(line);
                getline(readsFile, line);
                reads = boost::trim_copy(line);
                getline(readsFile, line);
                option = boost::trim_copy(line);
                getline(readsFile, line);
                score = boost::trim_copy(line);
            }

            for(int st = 0; st + readLength <= reads.length(); st += readLength) {
                nRealReads ++;

                bool filter = false;
                double qual = 0.0;
                for(int i = st; i < st + readLength; i ++) {
                    if(reads[i] == 'N' || reads[i] == 'n') {
                        filter == true;
                        break;
                    }
                    if(score[i] - '!' < BPQualBound) {
                        filter == true;
                        break;
                    }
                    qual += score[i] - '!';
                }

                if(filter || qual < readLength * AVGQualBound) {
                    nFilterOut ++;
                    continue;
                }

                if(pairId == 0) {
                    long long id = 0;
                    int p = 0;
                    for(; p < K - 1; p ++) {
                        char c = reads[st + p];
                        id = (id << 2) | parseBP(c);
                    }
                    for(; p < readLength; p ++) {
                        char c = reads[st + p];
                        id = ((id << 2) | parseBP(c)) & mask;

                        int pre = nFailKmer;
                        updateKmerCount(id, 1, nFailKmer);
                        //if(nFailKmer > pre) {
                        //    vector <string> rnaname;
                        //    boost::split(rnaname, identifier, boost::is_any_of(":"));
                        //    if(failset.count(rnaname[0]) > 0) {
                        //        failset[rnaname[0]] += 1;
                        //    } else {
                        //        failset[rnaname[0]] = 1;
                        //    }
                        //    if(rnaname[0] == "@NM_001320261") {
                        //        cout << endl;
                        //        cout << identifier << endl;
                        //        cout << id << endl;
                        //        cout << reads << endl;
                        //        cout << reads.substr(st +p, K) << endl;
                        //    } 
                        //}
                    }
                } else if (pairId == 1){
                    //cout << "Should not come here!" << endl;
                    long long antiId = 0;
                    int p = 0;
                    for(; p < K - 1; p ++) {
                        char c = reads[st + p];
                        antiId = (antiId >> 2) | (parseAntiBP(c) << shiftLeftMost);
                    }
                    for(; p < readLength; p ++) {
                        char c = reads[st + p];
                        antiId = ((antiId >> 2) | (parseAntiBP(c) << shiftLeftMost)) & mask;
                        updateKmerCount(antiId, 1, nFailKmer);
                    }
                } else {
                    cout << "Error: The number of reads file will never larger than 2!" << endl;
                    exit(-1);
                }
            }

            nRawReads ++;
            if(nRawReads % 1000000 == 0) {
                cout << "+++ " << nRawReads / 1000000 << ",000,000 reads loaded ..." << endl;
            }
        }

        readsFile.close();
    }

    //cout << failset.size() << endl;

    //for(auto i: failset) {
    //    cout << i.first << " " << i.second << endl;
    //}

    //exit(0);
    //getchar();

    cout << "=== Totally " << nRawReads << " raw reads loaded." << endl;
    cout << "=== Totally " << nRealReads << " real reads (with theoretical length) loaded." << endl;
    cout << "=== Reserve " << nRealReads - nFilterOut << " high quality real reads (with theoretical length). (Reserve " 
        << 100 - (double)nFilterOut / nRealReads * 100 << "%)" << endl;
    cout << "=== " << nFailKmer << " kmers failed to match with theoretical kmer profile. (Reserve " 
        << 100 - (double)nFailKmer / (nRealReads * (readLength - K + 1)) * 100 << "%)" << endl;

    cout << "=== Finally collect " << kmerTable.size() << " different possible kmers on Genome."<< endl;
    cout << "=== Finally collect " << kmerCount.size() << " different kmers occur in read."<< endl;
    cout << "=== About " << (nRealReads * (readLength - K + 1) - nFailKmer) / kmerCount.size() << " occurrences of kmer on average." << endl;
    cout << "*** Finish loading reads!" << endl;
}

void KmerHash::readGenome(string exonBoundaryPath, string predictExonPath) {
    vector <string> fileInfo;
    boost::split(fileInfo, exonBoundaryPath, boost::is_any_of("."));
    string ext = fileInfo.back();
    boost::algorithm::to_lower(ext);
    if(ext == "bed") {
        readFromBed(exonBoundaryPath, predictExonPath);
    } else {
        cout << "Error: Currently only support .BED format exon annotation file..." << endl;
        exit(-1);
    }
    return;
}

void KmerHash::readFromBed(string exonBoundaryPath, string predictExonPath) {
    const int COL_CHRNAME = 0;
    const int COL_GENEST = 1;
    const int COL_GENEED = 2;
    const int COL_STRAND = 5;
    const int COL_EXONNUM = 9;
    const int COL_EXONST = 11;
    const int COL_EXONLEN = 10;
    cout << "\n*** Start to load exon boundary annotation ..." << endl;

    ifstream exonBoundaryFile(exonBoundaryPath.c_str(), ios::in);
    string entry;
    string preChrName = "";
    int k = 0;
    vector <vector <int> > predictExon; 

    while(getline(exonBoundaryFile, entry)) {
        entry = boost::trim_copy(entry);
        vector <string> subcol;
        vector <string> subst;
        vector <string> sublen;
        boost::split(subcol, entry, boost::is_any_of("\t"));
        boost::split(subst, subcol[COL_EXONST], boost::is_any_of(","));
        boost::split(sublen, subcol[COL_EXONLEN], boost::is_any_of(","));

        vector <pair <int, int> > exons;
        vector <int> predictTemp;
        int exonnum = boost::lexical_cast<int> (subcol[COL_EXONNUM]);
        int genest = boost::lexical_cast<int> (subcol[COL_GENEST]);

        for(int i = 0; i < exonnum; i++) {
            int exonst = genest + boost::lexical_cast<int>(subst[i]);
            int exoned = exonst + boost::lexical_cast<int>(sublen[i]);

            if(exoned - exonst >= readLength) {
                exons.push_back(make_pair(exonst, exoned));
                predictTemp.push_back(i);
            } else {
                cout << "Caution: The length of exon " << i << " from gene " << geneBoundary.size() << " is shorter than reads." << endl;
            }
        }
        predictExon.push_back(predictTemp);

        geneBoundary.push_back(exons);
        string strand = subcol[COL_STRAND];
        geneStrand.push_back(strand);
        string chrName = subcol[COL_CHRNAME];
        if(chrName != preChrName) {
            preChrName = chrName;
            chromST.push_back(k);
            chromName.push_back(preChrName);
        }
        k ++;
    }
    exonBoundaryFile.close();

    NG = geneBoundary.size();
    for(int g = 0; g < NG; g++) {
        NE.push_back(geneBoundary[g].size());
    }
    cout << "=== There are " << NG << " genes from " << chromName.size() << " chromosomes" << endl;
    cout << "    (If the number of chromosomes is not right, please check whether the annotation is sorted.)" << endl;
    cout << "*** Finish loading exon boundary annotation !" << endl;
}

inline double KmerHash::positiveContribution(int st, int ed, int l, int r, int len) {
    int cil = min(len - readLength + 1, readLength - K + 1);
    int ret = min(l - st + 1, ed - r + 1);
    return min(ret, cil);
}

inline double KmerHash::negativeContribution(int st, int ed, int l, int r, int len) {
    int cil = min(readLength - 1 - len, readLength - K + 1);
    int ret = min(l - st + 1, ed - r + 1);
    return -min(ret, cil);
}

inline void KmerHash::insertKmerTable(long long id, long long exonid, double contr) {
    if(kmerTable.count(id) > 0) {
        if(kmerTable[id].count(exonid) > 0) {
            kmerTable[id][exonid] = kmerTable[id][exonid] + contr;
        } else {
            kmerTable[id][exonid] = contr;
        }
    } else {
        kmerTable[id] = unordered_map <long long, double> ();
        kmerTable[id][exonid] = contr;
    }
}

void KmerHash::buildKmerTable(string genomePath) {
    int chrId = 0;
    string chrName = "";
    string geneSeq = "";

    cout << "\n*** Start to build theoretical kmer profile ..." << endl;

    for(int g = 0; g < NG; g++) {
        if(chrId < chromST.size() && g >= chromST[chrId]) {
            chrName = chromName[chrId];
            geneSeq = "";
            cout << "+++ Build theoretical kmer profile from " << chrName << " ..." << endl;

            ifstream genomeFile;
            char buf[(1<<16)];
            genomeFile.rdbuf() -> pubsetbuf(buf, sizeof(buf));
            genomeFile.open((genomePath + "/" + chrName + ".fa").c_str(), ios::in);
            
            for(string line; getline(genomeFile, line);) {
                if(line[0] != '>') {
                    geneSeq += boost::trim_copy(line);
                }
            }
            genomeFile.close();
            chrId ++;
        }

        string strand = geneStrand[g];

        for(int e = 0; e < NE[g]; e++) {
            long long exonid;
            if(strand == "+") {
                exonid = encryptExonId(g, e, e);
            } else {
                exonid = encryptExonId(g, NE[g] - e - 1, NE[g] - e - 1);
            }

            int exonst = geneBoundary[g][e].first;
            int exoned = geneBoundary[g][e].second;
            int exonlen = exoned - exonst;

            if(exonlen >= readLength) {
                int p = 0;
                long long id = 0, antiId = 0;
                double tot = (exonlen - readLength + 1) * (readLength - K + 1);
                //double tot = (readLength - K + 1);

                for(; p < K - 1; p ++) {
                    char c = geneSeq[exonst + p];
                    id = (id << 2) | parseBP(c);
                    antiId = (antiId >> 2) | (parseAntiBP(c) << shiftLeftMost);
                }

                for(; exonst + p < exoned; p ++) {
                    char c = geneSeq[exonst + p];
                    id = ((id << 2) | parseBP(c)) & mask;
                    antiId = ((antiId >> 2) | (parseAntiBP(c) << shiftLeftMost)) & mask;

                    double contr = positiveContribution(0, exonlen, p - K + 1, p + 1, exonlen) / tot;

                    if(!strandSpec || strand == "+") {
                        insertKmerTable(id, exonid, contr);
                    }

                    if(!strandSpec || strand == "-") {
                        insertKmerTable(antiId, exonid, contr);
                    }
                }

            }
        }

        for(int ei = 0; ei < NE[g]; ei ++) {
            for(int ej = ei + 1; ej < NE[g]; ej ++) {
                long long exonid;
                if(strand == "+") {
                    exonid = encryptExonId(g, ei, ej);
                } else {
                    exonid = encryptExonId(g, NE[g] - ej - 1, NE[g] - ei - 1);
                }

                int exonedi = geneBoundary[g][ei].second;
                int exonstj = geneBoundary[g][ej].first;

                int p = 0;
                long long id = 0, antiId = 0;
                double tot = (readLength - 1) * (readLength - K + 1);
                //double tot = (readLength - K + 1);

                for(; p < K - 1; p ++) {
                    char c = geneSeq[exonedi - readLength + 1 + p];
                    id = (id << 2) | parseBP(c);
                    antiId = (antiId >> 2) | (parseAntiBP(c) << shiftLeftMost);
                }

                for(; p < readLength - 1; p ++) {
                    char c = geneSeq[exonedi - readLength + 1 + p];
                    id = ((id << 2) | parseBP(c)) & mask;
                    antiId = ((antiId >> 2) | (parseAntiBP(c) << shiftLeftMost)) & mask;

                    double contr = positiveContribution(0, 2*readLength-2, p-K+1, p+1, 2*readLength-2) / tot;

                    if(!strandSpec || strand == "+") {
                        insertKmerTable(id, exonid, contr);
                    }

                    if(!strandSpec || strand == "-") {
                        insertKmerTable(antiId, exonid, contr);
                    }
                }

                for(p = 0; p < readLength - 1; p ++) {
                    char c = geneSeq[exonstj + p];
                    id = ((id << 2) | parseBP(c)) & mask;
                    antiId = ((antiId >> 2) | (parseAntiBP(c) << shiftLeftMost)) & mask;

                    double contr = positiveContribution(0, 2*readLength-2, readLength+p-K, readLength+p, 2*readLength-2) / tot;

                    if(!strandSpec || strand == "+") {
                        insertKmerTable(id, exonid, contr);
                    }

                    if(!strandSpec || strand == "-") {
                        insertKmerTable(antiId, exonid, contr);
                    }
                }
            }
        }

    }
    cout << "*** Finish building theoretical kmer profile!" << endl;
}

template <typename Container> 
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};

void KmerHash::mergeKmerTable() {
    cout << "\n*** Start to merge kmers which share same profile ..." << endl;

    long long nEmptyTheoKmer = 0;
    long long nLowExpTheoKmer = 0;
    long long nSameKmers = 0;
    long long nUniqueKmers = 0;
    long long nTotal = 0;
    unordered_map <vector <long long>, long long,
                  container_hash<vector <long long> > > uniqueKmerTable;

    for(auto it = kmerTable.begin(); it != kmerTable.end(); ) {
        long long kmerId = it -> first;
        nTotal ++;
        if(nTotal % 10000000 == 0) {
            cout << "+++ " << nTotal / 1000000 << ",000,000 kmers processed ..." << endl;
        }

        if(kmerCount.count(kmerId) == 0) {
            it = kmerTable.erase(it);
            nEmptyTheoKmer ++;
            continue;
        }

        if(kmerCount[kmerId] <= lowExpKmerBound) {
            it = kmerTable.erase(it);
            kmerCount.erase(kmerId);
            nLowExpTheoKmer ++;
            continue;
        }

        vector <long long> exonIdList;
        for(auto exonProfile: it -> second) {
            exonIdList.push_back(exonProfile.first);
        }
        sort(exonIdList.begin(), exonIdList.end());

        if(uniqueKmerTable.count(exonIdList) > 0) {
            long long presentorKmerId = uniqueKmerTable[exonIdList];
            kmerCount[presentorKmerId] += kmerCount[kmerId];
            for(auto exon: exonIdList) {
                kmerTable[presentorKmerId][exon] += it -> second[exon];
            }
            kmerCount.erase(kmerId);
            it = kmerTable.erase(it);
            nSameKmers ++;
        } else {
            uniqueKmerTable[exonIdList] = kmerId;
            it ++;
            nUniqueKmers ++;
        }
    }

    NW = uniqueKmerTable.size();

    cout << "=== " << nEmptyTheoKmer << " kmers not occur in reads. ("
        << (double)nEmptyTheoKmer / nTotal * 100 << "%)." << endl;
    cout << "=== " << nLowExpTheoKmer << " kmers occur less than " << lowExpKmerBound 
        << " times in reads. (" << (double)nLowExpTheoKmer / nTotal * 100 << "%)." << endl;
    cout << "=== " << nSameKmers << " kmers share same profile. ("
        << (double)nSameKmers / nTotal * 100 << "%)." << endl;

    cout << "+++ Rehash kmer profile ..." << endl;
    double loadFactorKmerCount = kmerCount.load_factor();
    double loadFactorKmerTable = kmerTable.load_factor();
    kmerCount.rehash(kmerCount.size());
    kmerTable.rehash(kmerTable.size());
    cout << "    KmerCount load factor: Now " << kmerCount.load_factor() 
        << ", before " << loadFactorKmerCount << endl;
    cout << "    KmerTable load factor: Now " << kmerTable.load_factor() 
        << ", before " << loadFactorKmerTable << endl;
    cout << "=== " << NW << " kmers reserve after merging! " << endl;
    cout << "    KmerCountSize: " << kmerCount.size() << endl;
    cout << "    KmerTableSize: " << kmerTable.size() << endl;
    cout << "    (Shrinkage rate is " << (double)nUniqueKmers / nTotal * 100 << "%)." << endl;

    long long nExonUniqueKmers = 0;
    for(auto it: kmerTable) {
        if(it.second.size() == 1) {
            nExonUniqueKmers ++;
        }
    }
    cout << "=== " << nExonUniqueKmers << " kmers are originated from unique exon. (" 
        << (double) nExonUniqueKmers / NW * 100 << "%)." << endl;

    cout << "*** Finish merging kmer profile !" << endl;
}
