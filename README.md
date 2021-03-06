# Opal

*This work has been supported in part by Croatian Science Fundation under the project UIP-11-2013-7353.*

Opal (ex Swimd) is SIMD C/C++ library for massive optimal sequence alignment.
Opal is implemented mainly by Rognes's "Faster Smith-Waterman database searches with inter-sequence SIMD parallelisation". 
Main difference is that Opal offers support for AVX2 and 4 alignment modes instead of just Smith-Waterman.

#### Requirements
SSE4.1 or higher.
If AVX2 is available, Opal will consume two times more sequences and will therefore work two times faster.
By compiling code with makefile and running ./test, you can see if you have SSE4.1 or AVX2 and also see if everything is working.

#### Alignment modes
Opal offers 4 different modes of alignment: 1 local and 3 global modes, explained below.
* SW: Local alignment (Smith-Waterman). Useful for finding similar regions in not necessarily similar sequences, and also for finding shorter sequence in longer sequence (text searching).
* NW: Global alignment (Needleman-Wunsch). Useful for detecting if two sequences of similar lengths are similar.
* HW: Semi-global alignment. Insertions before query start and insertions after query end are not penalized. Useful for finding query in target.
* OV: Semi-global alignment. Insertions before start of either query or target and insertions after end of either query or target are not penalized. Useful when sequences are partially overlaping or one of sequences is contained in another.

#### Usage
To use Opal you just have to include opal.h in your code and compile your code together with opal.cpp using appropriate compiler flag for SSE, for example -msse4.1 (or just use -march which will detect your arhitecture and will also use appropriate SSE flag).  
Opal is written for C++11 standard, so you should make sure that you compile it according to that. For `gcc`, add flag `-std=c++11`.

```
...
#include "opal.h"
...
```

```
...
// Basic parameters
int alphabetLength = 4;
int gapOpen = 3;
int gapExt = 1;
int matchExt = 0;
int scoreMatrix[16] = {
    2, -1, -3, 0,
    -1, 4, -5, -1,
    -3, -5, 1, -10,
    0, -1, -10, 4
};

// Query
int queryLength = 10;
unsigned char query[10] = {0,1,3,2,1,0,3,0,1,1};

// Database
int dbLength = 4;
unsigned char dbSeq1[14] = {1,3,2,3,0,0,1,0,2,2,1,2,3,2};
unsigned char dbSeq2[12] = {2,1,1,3,2,0,0,2,2,0,2,1};
unsigned char dbSeq3[13] = {0,0,2,1,0,3,1,1,2,3,2,1,0};
unsigned char dbSeq4[9] = {2,3,3,3,1,1,2,2,0};
unsigned char* db[4] = { dbSeq1, dbSeq2, dbSeq3, dbSeq4 };
int dbSeqsLengths[4] = {14, 12, 13, 9};

// Results for each sequence in database
OpalSearchResult results[4];
for (int i = 0; i < 4; i++) {
    opalInitSearchResult(results[i]);
}

// Do calculation!
int resultCode = opalSearchDatabase(query, queryLength, db, dbLength, dbSeqsLengths,
                                    gapOpen, gapExt, matchExt, scoreMatrix, alphabetLength, results,
                                    OPAL_MODE_SW, OPAL_OVERFLOW_BUCKETS);

// Print scores
printf("%d %d %d %d\n", results[0].score, results[1].score, results[2].score, results[3].score);
...
```

For more examples of usage take a look at **test.cpp** and **opal_aligner.cpp**.
For detailed documentation check out **opal.h**.

## Opal aligner
In order to compile and use simple aligner that uses Opal run makefile in src:

    cd src
    make

Type `./opal_aligner` for help.

Examples of usage:

    ./opal_aligner ../test_data/query/O74807.fasta ../test_data/db/uniprot_sprot15.fasta
    ./opal_aligner -p ../test_data/query/O74807.fasta ../test_data/db/uniprot_sprot15.fasta
    ./opal_aligner -s ../test_data/query/P18080.fasta ../test_data/db/uniprot_sprot12071.fasta
    ./opal_aligner -s -a NW ../test_data/query/P18080.fasta ../test_data/db/uniprot_sprot12071.fasta

#### Test data
In test_data/ there are two directories, query/ and db/.
In query/ are fasta files with one sequence each, while in db/ are fasta files with many sequences each.
All fasta files were obtained from www.uniprot.org.
File in db/ with name uniprot_sprotN.fasta was generated by taking first N sequences from uniprot_sprot.fasta, which can be found at www.uniprot.org/downloads -> UniProtKB/Swiss-Prot.
