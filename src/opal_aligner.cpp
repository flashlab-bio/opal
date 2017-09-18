#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <sstream>
#include <climits>
#include <numeric>
#include <algorithm>
#include "opal.h"
#include "ScoreMatrix.hpp"

using namespace std;

#define DER   "\033[39m"        /* Reset color */
#define RED   "\033[31m"        /* Red */
#define UND   "\033[4m"         /* Underline */
#define DNU   "\033[24m"        /* Reset Underline */


bool readFastaSequences(unsigned char* seqin, unsigned char* alphabet, int alphabetLength, vector< vector<unsigned char> >* seqs, vector<string>* ids, bool rc);
bool readFastaSequences(FILE* &file, unsigned char* alphabet, int alphabetLength, vector< vector<unsigned char> >* seqs, vector<string>* ids, bool rc);
void printAlignment(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const OpalSearchResult result, const unsigned char* alphabet, const int wrap);

int main(int argc, char * const argv[]) {
    int matchScore = 0;
    int mismatchScore = 0;
    int gapOpen = 3;
    int gapExt = 1;
    int maxRes = 2;
    int wrap = 50;
    int matchExt = 0;
    ScoreMatrix scoreMatrix;

    //----------------------------- PARSE COMMAND LINE ------------------------//
    char matrixPattern[512];
    bool logging = false;
    bool revcomp = false;
    bool ispath = false;
    int modeCode = OPAL_MODE_SW;
    int searchType = OPAL_SEARCH_ALIGNMENT;
    int option;
    int matrixType = 0;
    while ((option = getopt(argc, argv, "a:o:e:b:n:w:m:x:hrp")) >= 0) {
        switch (option) {
        case 'a': modeCode = atoi(optarg); break;
        case 'o': gapOpen = atoi(optarg); break;
        case 'e': gapExt = atoi(optarg); break;
        case 'b': matchExt = atoi(optarg); break;
        case 'n': maxRes = atoi(optarg); break;
        case 'w': wrap = atoi(optarg); break;
        case 'm': switch (optarg[0]){
                    case '@': matrixType=1; break;
                    case ':': matrixType=2; break;
                    case 'n': matrixType=3; break;
                    case 'p': matrixType=4; break;
                  }
                  strcpy(matrixPattern, optarg+1); break;
        case 'x': searchType = atoi(optarg); break;
        case 'h': logging = true; break;
        case 'r': revcomp = true; break;
        case 'p': ispath = true; break;
        }
    }
    if (optind + 2 != argc) {
        if (!logging) return 1;
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: opal_aligner [options...] <query> <db>\n");
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -o NN\tis gap opening penalty. [default: 3]\n");
        fprintf(stderr, "  -e N\tN is gap extension penalty. [default: 1]\n"
                        "     Gap of length n will have penalty of g + (n - 1) * e.\n");
        fprintf(stderr, "  -n N\tN is max number of result entries.[default: 2]\n");
        fprintf(stderr, "  -w N\tN is wrap number of alignment view.[default: 50]\n");
        fprintf(stderr, "  -h\tIf set, more info will be output (logging mode).\n");
        fprintf(stderr, "  -r\tIf set, consider reverse-complement of query(s)\n");
        fprintf(stderr, "  -p\tIf set, treat input as path instead of sequences\n");
        fprintf(stderr, "  -a N\tAlignment mode that will be used,SW(3)|NW(0)|HW(1)|OV(2). [default: 3]\n");
        fprintf(stderr, "  -b N\tN is match extension bonus. [default: 0]\n"
                        "    This bonus is awarded to match following another match.\n");
        fprintf(stderr,
                "  -x\tsearch_level  Following search levels are available [default: %d]:\n"
                "    \t              %d - score\n"
                "    \t              %d - score, end location\n"
                "    \t              %d - score, end and start location and alignment\n",
                OPAL_SEARCH_ALIGNMENT, OPAL_SEARCH_SCORE, OPAL_SEARCH_SCORE_END, OPAL_SEARCH_ALIGNMENT);
        fprintf(stderr,
                "  -m\tpattern       Score matrix pattern. [default: n]:\n"
                "    \t              n       - nucleotide\n"
                "    \t              p       - protein\n"
                "    \t              @{path} - matrix filepath\n"
                "    \t              :{N[N]} - match score and mismatch score(optional)\n");
        return 1;
    }
    //-------------------------------------------------------------------------//

    // Set score matrix
    switch (matrixType){
        case 1: scoreMatrix = ScoreMatrix(matrixPattern); break;
        case 2: matchScore = matrixPattern[0] - 48;
                mismatchScore = matchScore;
                if (strlen(matrixPattern)>1) mismatchScore = matrixPattern[1] - 48;
                scoreMatrix = ScoreMatrix(matchScore, mismatchScore); break;
        case 3: case 0: scoreMatrix = ScoreMatrix::getBase50(); break;
        case 4: scoreMatrix = ScoreMatrix::getBlosum50(); break;
    }

    unsigned char* alphabet = scoreMatrix.getAlphabet();
    int alphabetLength = scoreMatrix.getAlphabetLength();
    int matrixLength = scoreMatrix.getMatrixLength();
    if (alphabetLength == 0 || matrixLength != alphabetLength * alphabetLength) {
        fprintf(stderr, "Invalid matrix!\n");
        return 1;
    }

    // Detect mode
    if (modeCode<0 || modeCode>3) {
        fprintf(stderr, "Invalid mode!\n");
        return 1;
    }
    if (logging) printf("Using (%d) alignment mode.\n", modeCode);

    // Build query
    vector< vector<unsigned char> >* querySequences = new vector< vector<unsigned char> >();
    vector<string>* queryID = new vector<string>();
    if (logging) printf("Reading query fasta file...\n");
    FILE* queryFile;
    if (ispath) {
        char* queryFilepath = argv[optind];
        queryFile = fopen(queryFilepath, "r");
        if (queryFile == 0) {
            fprintf(stderr, "Error: There is no file with name %s\n", queryFilepath);
            return 1;
        }
        readFastaSequences(queryFile, alphabet, alphabetLength, querySequences, queryID, revcomp);
    }else{
        readFastaSequences(reinterpret_cast<unsigned char *>(argv[optind]), alphabet, alphabetLength, querySequences, queryID, revcomp);
    }

    int qLength = queryID->size();
    if (qLength > 1) { maxRes = 1; }
    unsigned char* query;
    int queryLength = 0;
    if (ispath) { fclose(queryFile); }


    // Build db
    FILE* dbFile;
    if (ispath) {
        char* dbFilepath = argv[optind+1];
        dbFile = fopen(dbFilepath, "r");
        if (dbFile == 0) {
            fprintf(stderr, "Error: There is no file with name %s\n", dbFilepath);
            return 1;
        }
    }
    double cpuTime = 0;
    bool wholeDbRead = false;
    int dbTotalNumResidues = 0;  // Sum of lengths of all database sequences.
    int dbTotalLength = 0;  // Number of sequences in the database.
    vector<string>* dbID = new vector<string>(); // 
    while (!wholeDbRead) {
        vector< vector<unsigned char> >* dbSequences = new vector< vector<unsigned char> >();
        if (logging) printf("\nReading database fasta file...\n");
        // Chunk of database is read and processed (if database is not huge, there will be only one chunk).
        // We do this because if database is huge, it may not fit into memory.
        if (ispath) {
            wholeDbRead = readFastaSequences(dbFile, alphabet, alphabetLength, dbSequences, dbID, false);
        }else{
            wholeDbRead = readFastaSequences(reinterpret_cast<unsigned char *>(argv[optind+1]), alphabet, alphabetLength, dbSequences, dbID, false);
        }
        int dbLength = dbID->size();
        unsigned char** db = new unsigned char*[dbLength];
        int* dbSeqLengths = new int[dbLength];
        int dbNumResidues = 0;
        for (unsigned int i = 0; i < dbSequences->size(); i++) {
            db[i] = (*dbSequences)[i].data();
            dbSeqLengths[i] = (*dbSequences)[i].size();
            dbNumResidues += dbSeqLengths[i];
        }
        if (logging) printf("Read %d database sequences, %d residues total.\n", dbLength, dbNumResidues);

        dbTotalNumResidues += dbNumResidues;
        dbTotalLength += dbLength;
        if (wholeDbRead) {
            if (logging) printf("Whole database read: %d database sequences, %d residues in total.\n",
                   dbTotalLength, dbTotalNumResidues);
        }

        // ----------------------------- MAIN CALCULATION ----------------------------- //
        clock_t start = clock();
        OpalSearchResult** results = new OpalSearchResult*[dbLength];
        int resultCode;
        stringstream dumpRes;
        vector<int> idx(dbLength);
        iota(begin(idx), end(idx), 0);
        for (int j = 0; j < qLength; ++j) {
            queryLength = (*querySequences)[j].size();
            query = (*querySequences)[j].data();
            for (int i = 0; i < dbLength; i++) {
                results[i] = new OpalSearchResult;
                opalInitSearchResult(results[i]);
            }
            if (logging) printf("\nComparing query to database...");
            //fflush(stdout);
            resultCode = opalSearchDatabase(query, queryLength, db, dbLength, dbSeqLengths,
                                                 gapOpen, gapExt, matchExt, scoreMatrix.getMatrix(), alphabetLength,
                                                 results, searchType, modeCode, OPAL_OVERFLOW_BUCKETS);
            if (resultCode) {
                fprintf(stderr, "\nDatabase search failed with error code: %d\n", resultCode);
                continue;
            }
            // Sorting results
            sort(begin(idx), end(idx), [&](int i1, int i2) { return results[i1]->score > results[i2]->score; });
            // ---------------------------------------------------------------------------- //

            if (logging)
                printf("\n#<i>: <score> (<query start>, <target start>) (<query end>, <target end>)\n");
            for (int i = 0; i < min(dbLength, maxRes); i++) {
                printf(UND "%d\n" DNU, max(i,j));
                printf("::query(%d nt): %s::target(%d nt): %s::score %d", queryLength, (*queryID)[j].c_str(), dbSeqLengths[idx[i]], (*dbID)[idx[dbTotalLength - dbLength + i]].c_str(), results[idx[i]]->score);
                if (!logging)
                    dumpRes << "\n" << (*queryID)[j] << "\t" << (*dbID)[idx[dbTotalLength - dbLength + i]]
                            << "\t" << results[idx[i]]->startLocationQuery << "\t" << results[idx[i]]->endLocationQuery
                            << "\t" << results[idx[i]]->startLocationTarget << "\t" << results[idx[i]]->endLocationTarget;
                if (results[idx[i]]->startLocationQuery >= 0) {
                    printf(" (%d, %d)", results[idx[i]]->startLocationQuery, results[idx[i]]->startLocationTarget);
                } else {
                    printf(" (?, ?)");
                }
                if (results[idx[i]]->endLocationQuery >= 0) {
                    printf(" (%d, %d)", results[idx[i]]->endLocationQuery, results[idx[i]]->endLocationTarget);
                } else {
                    printf(" (?, ?)");
                }
                printf("\n");
                if (results[idx[i]]->alignment) {
                    printAlignment(query, queryLength, db[idx[i]], dbSeqLengths[idx[i]], *results[idx[i]], alphabet, wrap);
                }
            }
        }
        if (!logging) printf("%s", dumpRes.str().c_str());
        for (int i = 0; i < dbLength; i++) {
            if (results[i]->alignment) {
                free(results[i]->alignment);
            }
            delete (results[i]);
        }
        delete[] results;
        delete[] db;
        delete[] dbSeqLengths;
        delete dbSequences;

        clock_t finish = clock();
        cpuTime += ((double)(finish-start))/CLOCKS_PER_SEC;
        if (logging) printf("\nFinished!\n");
    }

    if (logging) printf("\nCpu time of searching: %.5lf\n", cpuTime);
    if (searchType != OPAL_SEARCH_ALIGNMENT && logging) {
        printf("GCUPS (giga cell updates per second): %.5lf\n",
               dbTotalNumResidues / 1000000000.0 * queryLength / cpuTime);
    }

    // Print this statistics only for SW because they are not valid for other modes.
    /*   if (!(strcmp(mode, "SW"))) {
        int for8, for16, for32;
        for8 = for16 = for32 = 0;
        double averageScore = 0;
        for (int i = 0; i < dbLength; i++) {
            averageScore += (double)scores[i] / dbLength;
            if (scores[i] < CHAR_MAX)
                for8++;
            else if (scores[i] < SHRT_MAX)
                for16++;
            else
                for32++;
        }
        printf("\nDatabase statistics:\n");
        printf("\tFor 8  (< %10d): %8d\n", CHAR_MAX, for8);
        printf("\tFor 16 (< %10d): %8d\n", SHRT_MAX, for16);
        printf("\tFor 32 (< %10d): %8d\n",  INT_MAX, for32);
        printf("\tAverage score: %lf\n", averageScore);
        }*/

    if (ispath) { fclose(dbFile); }
    // Free allocated space
    delete querySequences;
    fflush(stdout);
    return 0;
}




/** Reads sequences from fasta file. If it reads more than some amount of sequences, it will stop.
 * @param [in] file File pointer to database. It may not be the beginning of database.
 * @param [in] alphabet
 * @param [in] alphabetLength
 * @param [in] rc
 * @param [out] seqs Sequences will be stored here, each sequence as vector of indexes from alphabet.
 * @return true if reached end of file, otherwise false.
 */
bool readFastaSequences(FILE* &file, unsigned char* alphabet, int alphabetLength, vector< vector<unsigned char> >* seqs, vector<string>* ids, bool rc) {
    seqs->clear();

    unsigned char letterIdx[128]; //!< letterIdx[c] is index of letter c in alphabet
    for (int i = 0; i < alphabetLength; i++)
        if (alphabet[i] == '*') { // '*' represents all characters not in alphabet
            for (int j = 0; j < 128; j++)
                letterIdx[j] = i;
            break;
        }
    for (int i = 0; i < alphabetLength; i++) {
        letterIdx[alphabet[i]] = i;
        letterIdx[tolower(alphabet[i])] = i;
    }
    long numResiduesRead = 0;
    bool inHeader = false;
    bool inSequence = false;
    int buffSize = 4096;
    unsigned char buffer[buffSize];
    string buffID;
    int tmpPos;
    for(;;) {
        int read = fread(buffer, sizeof(char), buffSize, file);
        for (int i = 0; i < read; ++i) {
            unsigned char c = buffer[i];
            if (inHeader) { // I do nothing if in header
                if (c == '\n') {
                    inHeader = false;
                    ids->push_back(buffID);
                    buffID.clear();
                } else {
                    buffID += c;
                }
            } else {
                if (c == '>') {
                    if (seqs->size() > 0) {
                        numResiduesRead += seqs->back().size();
                        if (rc) {
                            seqs->push_back(vector<unsigned char>());
                            ids->push_back(ids->back());
                            ids->back().push_back(42);
                            tmpPos = seqs->size()-2;
                            for (int i=(*seqs)[tmpPos].size()-1; i > -1; i--) {
                                switch((*seqs)[tmpPos][i])
                                {
                                case 0:
                                    seqs->back().push_back(1);break;
                                case 1:
                                    seqs->back().push_back(0);break;
                                case 2:
                                    seqs->back().push_back(3);break;
                                case 3:
                                    seqs->back().push_back(2);break;
                                default:
                                    seqs->back().push_back((*seqs)[tmpPos][i]);
                                }
                            }
                        }
                        if (numResiduesRead > 1073741824L) {
                            fseek(file, i - read - 1, SEEK_CUR);
                            return false;
                        }
                    }
                    inHeader = true;
                    inSequence = false;
                } else {
                    if (c == '\r' || c == '\n')
                        continue;
                    // If starting new sequence, initialize it.
                    // Before that, check if we read more than 1GB of sequences,
                    // and if that is true, finish reading.
                    if (inSequence == false) {
                        inSequence = true;
                        seqs->push_back(vector<unsigned char>());
                    }
                    seqs->back().push_back(letterIdx[c]);
                }
            }
        }
        if (feof(file)) {
            if (ids->empty()) {
                ids->push_back("solo");
            }
            if (rc) {
                seqs->push_back(vector<unsigned char>());
                ids->push_back(ids->back());
                ids->back().push_back(42);
                for (int i=(*seqs)[seqs->size()-2].size()-1; i > -1; i--) {
                    switch((*seqs)[seqs->size()-2][i])
                    {
                        case 0:
                            seqs->back().push_back(1);break;
                        case 1:
                            seqs->back().push_back(0);break;
                        case 2:
                            seqs->back().push_back(3);break;
                        case 3:
                            seqs->back().push_back(2);break;
                    default:
                        seqs->back().push_back((*seqs)[seqs->size()-2][i]);
                    }
                }
            }
            break;
        }
    }
    return true;
}

bool readFastaSequences(unsigned char* seqin, unsigned char* alphabet, int alphabetLength, vector< vector<unsigned char> >* seqs, vector<string>* ids, bool rc) {
    seqs->clear();
    unsigned char letterIdx[128]; //!< letterIdx[c] is index of letter c in alphabet
    for (int i = 0; i < alphabetLength; i++)
        if (alphabet[i] == '*') { // '*' represents all characters not in alphabet
            for (int j = 0; j < 128; j++)
                letterIdx[j] = i;
            break;
        }
    for (int i = 0; i < alphabetLength; i++) {
        letterIdx[alphabet[i]] = i;
        letterIdx[tolower(alphabet[i])] = i;
    }
    bool inHeader = false;
    bool inSequence = false;
    string buffID;
    int tmpPos;
    while (*seqin) {
        if (inHeader) { // I do nothing if in header
            if (*seqin == '\n') {
                inHeader = false;
                ids->push_back(buffID);
                buffID.clear();
            } else {
                buffID += *seqin;
            }
        } else {
            if (*seqin == '>') {
                if (seqs->size() > 0) {
                    if (rc) {
                        seqs->push_back(vector<unsigned char>());
                        ids->push_back(ids->back());
                        ids->back().push_back(42);
                        tmpPos = seqs->size()-2;
                        for (int i=(*seqs)[tmpPos].size()-1; i > -1; i--) {
                            switch((*seqs)[tmpPos][i])
                            {
                            case 0:
                                seqs->back().push_back(1);break;
                            case 1:
                                seqs->back().push_back(0);break;
                            case 2:
                                seqs->back().push_back(3);break;
                            case 3:
                                seqs->back().push_back(2);break;
                            default:
                                seqs->back().push_back((*seqs)[tmpPos][i]);
                            }
                        }
                    }
                }
                inHeader = true;
                inSequence = false;
            } else {
                if (*seqin == '\r' || *seqin == '\n') {
                    ++seqin;
                    continue;
                }
                // If starting new sequence, initialize it.
                // Before that, check if we read more than 1GB of sequences,
                // and if that is true, finish reading.
                if (inSequence == false) {
                    inSequence = true;
                    seqs->push_back(vector<unsigned char>());
                }
                seqs->back().push_back(letterIdx[*seqin]);
            }
        }
        ++seqin;
    }
    if (ids->empty()) {
        ids->push_back("solo");
    }
    if (rc) {
        seqs->push_back(vector<unsigned char>());
        ids->push_back(ids->back());
        ids->back().push_back(42);
        for (int i=(*seqs)[seqs->size()-2].size()-1; i > -1; i--) {
            switch((*seqs)[seqs->size()-2][i])
            {
                case 0:
                    seqs->back().push_back(1);break;
                case 1:
                    seqs->back().push_back(0);break;
                case 2:
                    seqs->back().push_back(3);break;
                case 3:
                    seqs->back().push_back(2);break;
            default:
                seqs->back().push_back((*seqs)[seqs->size()-2][i]);
            }
        }
    }
    return true;
}

void printAlignment(const unsigned char* query, const int queryLength,
                    const unsigned char* target, const int targetLength,
                    const OpalSearchResult result, const unsigned char* alphabet, const int wrap) {
    int tIdx = result.startLocationTarget;
    int qIdx = result.startLocationQuery;
    /* What is this for?
      if (modeCode == EDLIB_MODE_HW) {
        tIdx = position;
        for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] != OPAL_ALIGN_DEL)
                tIdx--;
        }
      }
    */
    for (int start = 0; start < result.alignmentLength; start += wrap) {
        // target
        printf("T: ");
        int startTIdx = tIdx;
        for (int j = start; j < start + wrap && j < result.alignmentLength; j++) {
            if (result.alignment[j] == OPAL_ALIGN_DEL)
                printf("-");
            else if (result.alignment[j] == OPAL_ALIGN_MISMATCH)
                printf(RED "%c" DER, alphabet[target[tIdx++]]);
            else
                printf("%c", alphabet[target[tIdx++]]);
        }
        printf(" (%d - %d)\n", max(startTIdx, 0), tIdx - 1);
        // query
        printf("Q: ");
        int startQIdx = qIdx;
        for (int j = start; j < start + wrap && j < result.alignmentLength; j++) {
            if (result.alignment[j] == OPAL_ALIGN_INS)
                printf("-");
            else if (result.alignment[j] == OPAL_ALIGN_MISMATCH)
                printf(RED "%c" DER, alphabet[query[qIdx++]]);
            else
                printf("%c", alphabet[query[qIdx++]]);
        }
        printf(" (%d - %d)\n", max(startQIdx, 0), qIdx - 1);
    }
}
