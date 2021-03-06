#ifndef SCORE_MATRIX_HPP
#define SCORE_MATRIX_HPP

#include <vector>

using namespace std;

class ScoreMatrix {
 private:
    vector<unsigned char> alphabet; //!< letters in same order as columns/rows in matrix
    vector<int> matrix; //!< All rows of matrix concatenated. Has length: alphabetLength * alphabetLength

 public:
    ScoreMatrix();
    ScoreMatrix(vector<unsigned char> alphabet, vector<int> matrix);
    /**
     * Format of matrix file:
     *  - first line contains letters from alphabet naming columns separated with spaces.
     *  - each next line is one row of matrix (integers separated with spaces).
     */
    ScoreMatrix(const char* filepath);
    ScoreMatrix(const int match_score, const int mismatch_penalty);
    int getAlphabetLength();
    unsigned char* getAlphabet();
    int getMatrixLength();
    int* getMatrix();

    static ScoreMatrix getBlosum50();
    static ScoreMatrix getBase50();

 private:
    static vector<unsigned char> getBlosumAlphabet();
    static vector<unsigned char> getBaseAlphabet();
};

#endif // SCORE_MATRIX_HPP
