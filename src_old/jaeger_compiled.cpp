// [[Rcpp::plugins(openmp)]]  // Enable OpenMP in Rcpp
#include "edlib.h"
#include "RcppInt64.h"
#include <Rcpp.h>
#include <omp.h>  // Include OpenMP header for multi-threading
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdint.h>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <regex>
#include <unordered_map>
#include <cstdint>
#include <unordered_set>
#include <numeric>
#include <cmath>

using namespace std;
using namespace Rcpp;

typedef uint64_t Word;
static const int WORD_SIZE = sizeof(Word) * 8; // Size of Word in bits
static const Word WORD_1 = static_cast<Word>(1);
static const Word HIGH_BIT_MASK = WORD_1 << (WORD_SIZE - 1);  // 100..00
static const int MAX_UCHAR = 255;

// Data needed to find alignment.
struct AlignmentData {
  Word* Ps;
  Word* Ms;
  int* scores;
  int* firstBlocks;
  int* lastBlocks;
  
  AlignmentData(int maxNumBlocks, int targetLength) {
    // We build a complete table and mark first and last block for each column
    // (because algorithm is banded so only part of each columns is used).
    // TODO: do not build a whole table, but just enough blocks for each column.
    Ps     = new Word[maxNumBlocks * targetLength];
    Ms     = new Word[maxNumBlocks * targetLength];
    scores = new  int[maxNumBlocks * targetLength];
    firstBlocks = new int[targetLength];
    lastBlocks  = new int[targetLength];
  }
  
  ~AlignmentData() {
    delete[] Ps;
    delete[] Ms;
    delete[] scores;
    delete[] firstBlocks;
    delete[] lastBlocks;
  }
};

struct Block {
  Word P;  // Pvin
  Word M;  // Mvin
  int score; // score of last cell in block;
  
  Block() {}
  Block(Word p, Word m, int s) :P(p), M(m), score(s) {}
};


/**
 * Defines equality relation on alphabet characters.
 * By default each character is always equal only to itself, but you can also provide additional equalities.
 */
class EqualityDefinition {
private:
  bool matrix[MAX_UCHAR + 1][MAX_UCHAR + 1];
public:
  EqualityDefinition(const string& alphabet,
    const EdlibEqualityPair* additionalEqualities = NULL,
    const int additionalEqualitiesLength = 0) {
    for (int i = 0; i < static_cast<int>(alphabet.size()); i++) {
      for (int j = 0; j < static_cast<int>(alphabet.size()); j++) {
        matrix[i][j] = (i == j);
      }
    }
    if (additionalEqualities != NULL) {
      for (int i = 0; i < additionalEqualitiesLength; i++) {
        size_t firstTransformed = alphabet.find(additionalEqualities[i].first);
        size_t secondTransformed = alphabet.find(additionalEqualities[i].second);
        if (firstTransformed != string::npos && secondTransformed != string::npos) {
          matrix[firstTransformed][secondTransformed] = matrix[secondTransformed][firstTransformed] = true;
        }
      }
    }
  }
  
  /**
   * @param a  Element from transformed sequence.
   * @param b  Element from transformed sequence.
   * @return True if a and b are defined as equal, false otherwise.
   */
  bool areEqual(unsigned char a, unsigned char b) const {
    return matrix[a][b];
  }
};

static int myersCalcEditDistanceSemiGlobal(const Word* Peq, int W, int maxNumBlocks,
  int queryLength,
  const unsigned char* target, int targetLength,
  int k, EdlibAlignMode mode,
  int* bestScore_, int** positions_, int* numPositions_);

static int myersCalcEditDistanceNW(const Word* Peq, int W, int maxNumBlocks,
  int queryLength,
  const unsigned char* target, int targetLength,
  int k, int* bestScore_,
  int* position_, bool findAlignment,
  AlignmentData** alignData, int targetStopPosition);


static int obtainAlignment(
    const unsigned char* query, const unsigned char* rQuery, int queryLength,
    const unsigned char* target, const unsigned char* rTarget, int targetLength,
    const EqualityDefinition& equalityDefinition, int alphabetLength, int bestScore,
    unsigned char** alignment, int* alignmentLength);

static int obtainAlignmentHirschberg(
    const unsigned char* query, const unsigned char* rQuery, int queryLength,
    const unsigned char* target, const unsigned char* rTarget, int targetLength,
    const EqualityDefinition& equalityDefinition, int alphabetLength, int bestScore,
    unsigned char** alignment, int* alignmentLength);

static int obtainAlignmentTraceback(int queryLength, int targetLength,
  int bestScore, const AlignmentData* alignData,
  unsigned char** alignment, int* alignmentLength);

static string transformSequences(const char* queryOriginal, int queryLength,
  const char* targetOriginal, int targetLength,
  unsigned char** queryTransformed,
  unsigned char** targetTransformed);

static inline int ceilDiv(int x, int y);

static inline unsigned char* createReverseCopy(const unsigned char* seq, int length);

static inline Word* buildPeq(const int alphabetLength,
  const unsigned char* query,
  const int queryLength,
  const EqualityDefinition& equalityDefinition);


/**
 * Main edlib method.
 */
extern "C" EdlibAlignResult edlibAlign(const char* const queryOriginal, const int queryLength,
  const char* const targetOriginal, const int targetLength,
  const EdlibAlignConfig config) {
  EdlibAlignResult result;
  result.status = EDLIB_STATUS_OK;
  result.editDistance = -1;
  result.endLocations = result.startLocations = NULL;
  result.numLocations = 0;
  result.alignment = NULL;
  result.alignmentLength = 0;
  result.alphabetLength = 0;
  
  /*------------ TRANSFORM SEQUENCES AND RECOGNIZE ALPHABET -----------*/
  unsigned char* query, * target;
  string alphabet = transformSequences(queryOriginal, queryLength, targetOriginal, targetLength,
    &query, &target);
  result.alphabetLength = static_cast<int>(alphabet.size());
  /*-------------------------------------------------------*/
  
  // Handle special situation when at least one of the sequences has length 0.
  if (queryLength == 0 || targetLength == 0) {
    if (config.mode == EDLIB_MODE_NW) {
      result.editDistance = std::max(queryLength, targetLength);
      result.endLocations = static_cast<int *>(malloc(sizeof(int) * 1));
      result.endLocations[0] = targetLength - 1;
      result.numLocations = 1;
    } else if (config.mode == EDLIB_MODE_SHW || config.mode == EDLIB_MODE_HW) {
      result.editDistance = queryLength;
      result.endLocations = static_cast<int *>(malloc(sizeof(int) * 1));
      result.endLocations[0] = -1;
      result.numLocations = 1;
    } else {
      result.status = EDLIB_STATUS_ERROR;
    }
    
    free(query);
    free(target);
    return result;
  }
  
  /*--------------------- INITIALIZATION ------------------*/
  int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE); // bmax in Myers
  int W = maxNumBlocks * WORD_SIZE - queryLength; // number of redundant cells in last level blocks
  EqualityDefinition equalityDefinition(alphabet, config.additionalEqualities, config.additionalEqualitiesLength);
  Word* Peq = buildPeq(static_cast<int>(alphabet.size()), query, queryLength, equalityDefinition);
  /*-------------------------------------------------------*/
  
  /*------------------ MAIN CALCULATION -------------------*/
  // TODO: Store alignment data only after k is determined? That could make things faster.
  int positionNW; // Used only when mode is NW.
  AlignmentData* alignData = NULL;
  bool dynamicK = false;
  int k = config.k;
  if (k < 0) { // If valid k is not given, auto-adjust k until solution is found.
    dynamicK = true;
    k = WORD_SIZE; // Gives better results than smaller k.
  }
  
  do {
    if (config.mode == EDLIB_MODE_HW || config.mode == EDLIB_MODE_SHW) {
      myersCalcEditDistanceSemiGlobal(Peq, W, maxNumBlocks,
        queryLength, target, targetLength,
        k, config.mode, &(result.editDistance),
        &(result.endLocations), &(result.numLocations));
    } else {  // mode == EDLIB_MODE_NW
      myersCalcEditDistanceNW(Peq, W, maxNumBlocks,
        queryLength, target, targetLength,
        k, &(result.editDistance), &positionNW,
        false, &alignData, -1);
    }
    k *= 2;
  } while(dynamicK && result.editDistance == -1);
  
  if (result.editDistance >= 0) {  // If there is solution.
    // If NW mode, set end location explicitly.
    if (config.mode == EDLIB_MODE_NW) {
      result.endLocations = static_cast<int *>(malloc(sizeof(int) * 1));
      result.endLocations[0] = targetLength - 1;
      result.numLocations = 1;
    }
    
    // Find starting locations.
    if (config.task == EDLIB_TASK_LOC || config.task == EDLIB_TASK_PATH) {
      result.startLocations = static_cast<int *>(malloc(result.numLocations * sizeof(int)));
      if (config.mode == EDLIB_MODE_HW) {  // If HW, I need to calculate start locations.
        const unsigned char* rTarget = createReverseCopy(target, targetLength);
        const unsigned char* rQuery  = createReverseCopy(query, queryLength);
        // Peq for reversed query.
        Word* rPeq = buildPeq(static_cast<int>(alphabet.size()), rQuery, queryLength, equalityDefinition);
        for (int i = 0; i < result.numLocations; i++) {
          int endLocation = result.endLocations[i];
          if (endLocation == -1) {
            // NOTE: Sometimes one of optimal solutions is that query starts before target, like this:
            //                       AAGG <- target
            //                   CCTT     <- query
            //   It will never be only optimal solution and it does not happen often, however it is
            //   possible and in that case end location will be -1. What should we do with that?
            //   Should we just skip reporting such end location, although it is a solution?
            //   If we do report it, what is the start location? -4? -1? Nothing?
            // TODO: Figure this out. This has to do in general with how we think about start
            //   and end locations.
            //   Also, we have alignment later relying on this locations to limit the space of it's
            //   search -> how can it do it right if these locations are negative or incorrect?
            result.startLocations[i] = 0;  // I put 0 for now, but it does not make much sense.
          } else {
            int bestScoreSHW, numPositionsSHW;
            int* positionsSHW;
            myersCalcEditDistanceSemiGlobal(
              rPeq, W, maxNumBlocks,
              queryLength, rTarget + targetLength - endLocation - 1, endLocation + 1,
              result.editDistance, EDLIB_MODE_SHW,
              &bestScoreSHW, &positionsSHW, &numPositionsSHW);
            // Taking last location as start ensures that alignment will not start with insertions
            // if it can start with mismatches instead.
            result.startLocations[i] = endLocation - positionsSHW[numPositionsSHW - 1];
            free(positionsSHW);
          }
        }
        delete[] rTarget;
        delete[] rQuery;
        delete[] rPeq;
      } else {  // If mode is SHW or NW
        for (int i = 0; i < result.numLocations; i++) {
          result.startLocations[i] = 0;
        }
      }
    }
    
    // Find alignment -> all comes down to finding alignment for NW.
    // Currently we return alignment only for first pair of locations.
    if (config.task == EDLIB_TASK_PATH) {
      int alnStartLocation = result.startLocations[0];
      int alnEndLocation = result.endLocations[0];
      const unsigned char* alnTarget = target + alnStartLocation;
      const int alnTargetLength = alnEndLocation - alnStartLocation + 1;
      const unsigned char* rAlnTarget = createReverseCopy(alnTarget, alnTargetLength);
      const unsigned char* rQuery  = createReverseCopy(query, queryLength);
      obtainAlignment(query, rQuery, queryLength,
        alnTarget, rAlnTarget, alnTargetLength,
        equalityDefinition, static_cast<int>(alphabet.size()), result.editDistance,
        &(result.alignment), &(result.alignmentLength));
      delete[] rAlnTarget;
      delete[] rQuery;
    }
  }
  /*-------------------------------------------------------*/
  
  //--- Free memory ---//
  delete[] Peq;
  free(query);
  free(target);
  if (alignData) delete alignData;
  //-------------------//
  
  return result;
}

extern "C" char* edlibAlignmentToCigar(const unsigned char* const alignment, const int alignmentLength,
  const EdlibCigarFormat cigarFormat) {
  if (cigarFormat != EDLIB_CIGAR_EXTENDED && cigarFormat != EDLIB_CIGAR_STANDARD) {
    return 0;
  }
  
  // Maps move code from alignment to char in cigar.
  //                        0    1    2    3
  char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
  if (cigarFormat == EDLIB_CIGAR_STANDARD) {
    moveCodeToChar[0] = moveCodeToChar[3] = 'M';
  }
  
  vector<char>* cigar = new vector<char>();
  char lastMove = 0;  // Char of last move. 0 if there was no previous move.
  int numOfSameMoves = 0;
  for (int i = 0; i <= alignmentLength; i++) {
    // if new sequence of same moves started
    if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
      // Write number of moves to cigar string.
      int numDigits = 0;
      for (; numOfSameMoves; numOfSameMoves /= 10) {
        cigar->push_back('0' + numOfSameMoves % 10);
        numDigits++;
      }
      reverse(cigar->end() - numDigits, cigar->end());
      // Write code of move to cigar string.
      cigar->push_back(lastMove);
      // If not at the end, start new sequence of moves.
      if (i < alignmentLength) {
        // Check if alignment has valid values.
        if (alignment[i] > 3) {
          delete cigar;
          return 0;
        }
        numOfSameMoves = 0;
      }
    }
    if (i < alignmentLength) {
      lastMove = moveCodeToChar[alignment[i]];
      numOfSameMoves++;
    }
  }
  cigar->push_back(0);  // Null character termination.
  char* cigar_ = static_cast<char *>(malloc(cigar->size() * sizeof(char)));
  memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
  delete cigar;
  
  return cigar_;
}

/**
 * Build Peq table for given query and alphabet.
 * Peq is table of dimensions alphabetLength+1 x maxNumBlocks.
 * Bit i of Peq[s * maxNumBlocks + b] is 1 if i-th symbol from block b of query equals symbol s, otherwise it is 0.
 * NOTICE: free returned array with delete[]!
 */
static inline Word* buildPeq(const int alphabetLength,
  const unsigned char* const query,
  const int queryLength,
  const EqualityDefinition& equalityDefinition) {
  int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  // table of dimensions alphabetLength+1 x maxNumBlocks. Last symbol is wildcard.
  Word* Peq = new Word[(alphabetLength + 1) * maxNumBlocks];
  
  // Build Peq (1 is match, 0 is mismatch). NOTE: last column is wildcard(symbol that matches anything) with just 1s
  for (int symbol = 0; symbol <= alphabetLength; symbol++) {
    for (int b = 0; b < maxNumBlocks; b++) {
      if (symbol < alphabetLength) {
        Peq[symbol * maxNumBlocks + b] = 0;
        for (int r = (b+1) * WORD_SIZE - 1; r >= b * WORD_SIZE; r--) {
          Peq[symbol * maxNumBlocks + b] <<= 1;
          // NOTE: We pretend like query is padded at the end with W wildcard symbols
          if (r >= queryLength || equalityDefinition.areEqual(query[r], symbol))
            Peq[symbol * maxNumBlocks + b] += 1;
        }
      } else { // Last symbol is wildcard, so it is all 1s
        Peq[symbol * maxNumBlocks + b] = static_cast<Word>(-1);
      }
    }
  }
  
  return Peq;
}


/**
 * Returns new sequence that is reverse of given sequence.
 * Free returned array with delete[].
 */
static inline unsigned char* createReverseCopy(const unsigned char* const seq, const int length) {
  unsigned char* rSeq = new unsigned char[length];
  for (int i = 0; i < length; i++) {
    rSeq[i] = seq[length - i - 1];
  }
  return rSeq;
}

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word (one most to the left) is most bottom cell of block from column.
 * Pv[i] and Mv[i] define vin of cell[i]: vin = cell[i] - cell[i-1].
 * @param [in] Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
 * @param [in] Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.
 * @param [in] Eq  Bitset, Eq[i] == 1 if match, 0 if mismatch.
 * @param [in] hin  Will be +1, 0 or -1.
 * @param [out] PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
 * @param [out] MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
 * @param [out] hout  Will be +1, 0 or -1.
 */
static inline int calculateBlock(Word Pv, Word Mv, Word Eq, const int hin,
  Word &PvOut, Word &MvOut) {
  // hin can be 1, -1 or 0.
  // 1  -> 00...01
  // 0  -> 00...00
  // -1 -> 11...11 (2-complement)
  
  Word hinIsNeg = static_cast<Word>(hin >> 2) & WORD_1; // 00...001 if hin is -1, 00...000 if 0 or 1
  
  Word Xv = Eq | Mv;
  // This is instruction below written using 'if': if (hin < 0) Eq |= (Word)1;
  Eq |= hinIsNeg;
  Word Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
  
  Word Ph = Mv | ~(Xh | Pv);
  Word Mh = Pv & Xh;
  
  int hout = 0;
  // This is instruction below written using 'if': if (Ph & HIGH_BIT_MASK) hout = 1;
  hout = (Ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1);
  // This is instruction below written using 'if': if (Mh & HIGH_BIT_MASK) hout = -1;
  hout -= (Mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1);
  
  Ph <<= 1;
  Mh <<= 1;
  
  // This is instruction below written using 'if': if (hin < 0) Mh |= (Word)1;
  Mh |= hinIsNeg;
  // This is instruction below written using 'if': if (hin > 0) Ph |= (Word)1;
  Ph |= static_cast<Word>((hin + 1) >> 1);
  
  PvOut = Mh | ~(Xv | Ph);
  MvOut = Ph & Xv;
  
  return hout;
}

/**
 * Does ceiling division x / y.
 * Note: x and y must be non-negative and x + y must not overflow.
 */
static inline int ceilDiv(const int x, const int y) {
  return x % y ? x / y + 1 : x / y;
}

static inline int min(const int x, const int y) {
  return x < y ? x : y;
}

static inline int max(const int x, const int y) {
  return x > y ? x : y;
}


/**
 * @param [in] block
 * @return Values of cells in block, starting with bottom cell in block.
 */
static inline vector<int> getBlockCellValues(const Block block) {
  vector<int> scores(WORD_SIZE);
  int score = block.score;
  Word mask = HIGH_BIT_MASK;
  for (int i = 0; i < WORD_SIZE - 1; i++) {
    scores[i] = score;
    if (block.P & mask) score--;
    if (block.M & mask) score++;
    mask >>= 1;
  }
  scores[WORD_SIZE - 1] = score;
  return scores;
}

/**
 * Writes values of cells in block into given array, starting with first/top cell.
 * @param [in] block
 * @param [out] dest  Array into which cell values are written. Must have size of at least WORD_SIZE.
 */
static inline void readBlock(const Block block, int* const dest) {
  int score = block.score;
  Word mask = HIGH_BIT_MASK;
  for (int i = 0; i < WORD_SIZE - 1; i++) {
    dest[WORD_SIZE - 1 - i] = score;
    if (block.P & mask) score--;
    if (block.M & mask) score++;
    mask >>= 1;
  }
  dest[0] = score;
}

/**
 * Writes values of cells in block into given array, starting with last/bottom cell.
 * @param [in] block
 * @param [out] dest  Array into which cell values are written. Must have size of at least WORD_SIZE.
 */
static inline void readBlockReverse(const Block block, int* const dest) {
  int score = block.score;
  Word mask = HIGH_BIT_MASK;
  for (int i = 0; i < WORD_SIZE - 1; i++) {
    dest[i] = score;
    if (block.P & mask) score--;
    if (block.M & mask) score++;
    mask >>= 1;
  }
  dest[WORD_SIZE - 1] = score;
}

/**
 * @param [in] block
 * @param [in] k
 * @return True if all cells in block have value larger than k, otherwise false.
 */
static inline bool allBlockCellsLarger(const Block block, const int k) {
  vector<int> scores = getBlockCellValues(block);
  for (int i = 0; i < WORD_SIZE; i++) {
    if (scores[i] <= k) return false;
  }
  return true;
}


/**
 * Uses Myers' bit-vector algorithm to find edit distance for one of semi-global alignment methods.
 * @param [in] Peq  Query profile.
 * @param [in] W  Size of padding in last block.
 *                TODO: Calculate this directly from query, instead of passing it.
 * @param [in] maxNumBlocks  Number of blocks needed to cover the whole query.
 *                           TODO: Calculate this directly from query, instead of passing it.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] targetLength
 * @param [in] k
 * @param [in] mode  EDLIB_MODE_HW or EDLIB_MODE_SHW
 * @param [out] bestScore_  Edit distance.
 * @param [out] positions_  Array of 0-indexed positions in target at which best score was found.
 Make sure to free this array with free().
 * @param [out] numPositions_  Number of positions in the positions_ array.
 * @return Status.
 */
static int myersCalcEditDistanceSemiGlobal(
    const Word* const Peq, const int W, const int maxNumBlocks,
    const int queryLength,
    const unsigned char* const target, const int targetLength,
    int k, const EdlibAlignMode mode,
    int* const bestScore_, int** const positions_, int* const numPositions_) {
  *positions_ = NULL;
  *numPositions_ = 0;
  
  // firstBlock is 0-based index of first block in Ukkonen band.
  // lastBlock is 0-based index of last block in Ukkonen band.
  int firstBlock = 0;
  int lastBlock = min(ceilDiv(k + 1, WORD_SIZE), maxNumBlocks) - 1; // y in Myers
  Block *bl; // Current block
  
  Block* blocks = new Block[maxNumBlocks];
  
  // For HW, solution will never be larger then queryLength.
  if (mode == EDLIB_MODE_HW) {
    k = min(queryLength, k);
  }
  
  // Each STRONG_REDUCE_NUM column is reduced in more expensive way.
  // This gives speed up of about 2 times for small k.
  const int STRONG_REDUCE_NUM = 2048;
  
  // Initialize P, M and score
  bl = blocks;
  for (int b = 0; b <= lastBlock; b++) {
    bl->score = (b + 1) * WORD_SIZE;
    bl->P = static_cast<Word>(-1); // All 1s
    bl->M = static_cast<Word>(0);
    bl++;
  }
  
  int bestScore = -1;
  vector<int> positions; // TODO: Maybe put this on heap?
  const int startHout = mode == EDLIB_MODE_HW ? 0 : 1; // If 0 then gap before query is not penalized;
  const unsigned char* targetChar = target;
  for (int c = 0; c < targetLength; c++) { // for each column
    const Word* Peq_c = Peq + (*targetChar) * maxNumBlocks;
    
    //----------------------- Calculate column -------------------------//
    int hout = startHout;
    bl = blocks + firstBlock;
    Peq_c += firstBlock;
    for (int b = firstBlock; b <= lastBlock; b++) {
      hout = calculateBlock(bl->P, bl->M, *Peq_c, hout, bl->P, bl->M);
      bl->score += hout;
      bl++; Peq_c++;
    }
    bl--; Peq_c--;
    //------------------------------------------------------------------//
    
    //---------- Adjust number of blocks according to Ukkonen ----------//
    if ((lastBlock < maxNumBlocks - 1) && (bl->score - hout <= k) // bl is pointing to last block
      && ((*(Peq_c + 1) & WORD_1) || hout < 0)) { // Peq_c is pointing to last block
      // If score of left block is not too big, calculate one more block
      lastBlock++; bl++; Peq_c++;
      bl->P = static_cast<Word>(-1); // All 1s
      bl->M = static_cast<Word>(0);
      bl->score = (bl - 1)->score - hout + WORD_SIZE + calculateBlock(bl->P, bl->M, *Peq_c, hout, bl->P, bl->M);
    } else {
      while (lastBlock >= firstBlock && bl->score >= k + WORD_SIZE) {
        lastBlock--; bl--; Peq_c--;
      }
    }
    
    // Every some columns, do some expensive but also more efficient block reducing.
    // This is important!
    //
    // Reduce the band by decreasing last block if possible.
    if (c % STRONG_REDUCE_NUM == 0) {
      while (lastBlock >= 0 && lastBlock >= firstBlock && allBlockCellsLarger(*bl, k)) {
        lastBlock--; bl--; Peq_c--;
      }
    }
    // For HW, even if all cells are > k, there still may be solution in next
    // column because starting conditions at upper boundary are 0.
    // That means that first block is always candidate for solution,
    // and we can never end calculation before last column.
    if (mode == EDLIB_MODE_HW && lastBlock == -1) {
      lastBlock++; bl++; Peq_c++;
    }
    
    // Reduce band by increasing first block if possible. Not applicable to HW.
    if (mode != EDLIB_MODE_HW) {
      while (firstBlock <= lastBlock && blocks[firstBlock].score >= k + WORD_SIZE) {
        firstBlock++;
      }
      if (c % STRONG_REDUCE_NUM == 0) { // Do strong reduction every some blocks
        while (firstBlock <= lastBlock && allBlockCellsLarger(blocks[firstBlock], k)) {
          firstBlock++;
        }
      }
    }
    
    // If band stops to exist finish
    if (lastBlock < firstBlock) {
      *bestScore_ = bestScore;
      if (bestScore != -1) {
        *positions_ = static_cast<int *>(malloc(sizeof(int) * static_cast<int>(positions.size())));
        *numPositions_ = static_cast<int>(positions.size());
        copy(positions.begin(), positions.end(), *positions_);
      }
      delete[] blocks;
      return EDLIB_STATUS_OK;
    }
    //------------------------------------------------------------------//
    
    //------------------------- Update best score ----------------------//
    if (lastBlock == maxNumBlocks - 1) {
      int colScore = bl->score;
      if (colScore <= k) { // Scores > k dont have correct values (so we cannot use them), but are certainly > k.
        // NOTE: Score that I find in column c is actually score from column c-W
        if (bestScore == -1 || colScore <= bestScore) {
          if (colScore != bestScore) {
            positions.clear();
            bestScore = colScore;
            // Change k so we will look only for equal or better
            // scores then the best found so far.
            k = bestScore;
          }
          positions.push_back(c - W);
        }
      }
    }
    //------------------------------------------------------------------//
    
    targetChar++;
  }
  
  
  // Obtain results for last W columns from last column.
  if (lastBlock == maxNumBlocks - 1) {
    vector<int> blockScores = getBlockCellValues(*bl);
    for (int i = 0; i < W; i++) {
      int colScore = blockScores[i + 1];
      if (colScore <= k && (bestScore == -1 || colScore <= bestScore)) {
        if (colScore != bestScore) {
          positions.clear();
          k = bestScore = colScore;
        }
        positions.push_back(targetLength - W + i);
      }
    }
  }
  
  *bestScore_ = bestScore;
  if (bestScore != -1) {
    *positions_ = static_cast<int *>(malloc(sizeof(int) * static_cast<int>(positions.size())));
    *numPositions_ = static_cast<int>(positions.size());
    copy(positions.begin(), positions.end(), *positions_);
  }
  
  delete[] blocks;
  return EDLIB_STATUS_OK;
}


/**
 * Uses Myers' bit-vector algorithm to find edit distance for global(NW) alignment method.
 * @param [in] Peq  Query profile.
 * @param [in] W  Size of padding in last block.
 *                TODO: Calculate this directly from query, instead of passing it.
 * @param [in] maxNumBlocks  Number of blocks needed to cover the whole query.
 *                           TODO: Calculate this directly from query, instead of passing it.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] targetLength
 * @param [in] k
 * @param [out] bestScore_  Edit distance.
 * @param [out] position_  0-indexed position in target at which best score was found.
 * @param [in] findAlignment  If true, whole matrix is remembered and alignment data is returned.
 *                            Quadratic amount of memory is consumed.
 * @param [out] alignData  Data needed for alignment traceback (for reconstruction of alignment).
 *                         Set only if findAlignment is set to true, otherwise it is NULL.
 *                         Make sure to free this array using delete[].
 * @param [out] targetStopPosition  If set to -1, whole calculation is performed normally, as expected.
 *         If set to p, calculation is performed up to position p in target (inclusive)
 *         and column p is returned as the only column in alignData.
 * @return Status.
 */
static int myersCalcEditDistanceNW(const Word* const Peq, const int W, const int maxNumBlocks,
  const int queryLength,
  const unsigned char* const target, const int targetLength,
  int k, int* const bestScore_,
  int* const position_, const bool findAlignment,
  AlignmentData** const alignData, const int targetStopPosition) {
  if (targetStopPosition > -1 && findAlignment) {
    // They can not be both set at the same time!
    return EDLIB_STATUS_ERROR;
  }
  
  // Each STRONG_REDUCE_NUM column is reduced in more expensive way.
  const int STRONG_REDUCE_NUM = 2048; // TODO: Choose this number dinamically (based on query and target lengths?), so it does not affect speed of computation
  
  if (k < abs(targetLength - queryLength)) {
    *bestScore_ = *position_ = -1;
    return EDLIB_STATUS_OK;
  }
  
  k = min(k, max(queryLength, targetLength));  // Upper bound for k
  
  // firstBlock is 0-based index of first block in Ukkonen band.
  // lastBlock is 0-based index of last block in Ukkonen band.
  int firstBlock = 0;
  // This is optimal now, by my formula.
  int lastBlock = min(maxNumBlocks, ceilDiv(min(k, (k + queryLength - targetLength) / 2) + 1, WORD_SIZE)) - 1;
  Block* bl; // Current block
  
  Block* blocks = new Block[maxNumBlocks];
  
  // Initialize P, M and score
  bl = blocks;
  for (int b = 0; b <= lastBlock; b++) {
    bl->score = (b + 1) * WORD_SIZE;
    bl->P = static_cast<Word>(-1); // All 1s
    bl->M = static_cast<Word>(0);
    bl++;
  }
  
  // If we want to find alignment, we have to store needed data.
  if (findAlignment)
    *alignData = new AlignmentData(maxNumBlocks, targetLength);
  else if (targetStopPosition > -1)
    *alignData = new AlignmentData(maxNumBlocks, 1);
  else
    *alignData = NULL;
  
  const unsigned char* targetChar = target;
  for (int c = 0; c < targetLength; c++) { // for each column
    const Word* Peq_c = Peq + *targetChar * maxNumBlocks;
    
    //----------------------- Calculate column -------------------------//
    int hout = 1;
    bl = blocks + firstBlock;
    for (int b = firstBlock; b <= lastBlock; b++) {
      hout = calculateBlock(bl->P, bl->M, Peq_c[b], hout, bl->P, bl->M);
      bl->score += hout;
      bl++;
    }
    bl--;
    //------------------------------------------------------------------//
    // bl now points to last block
    
    // Update k. I do it only on end of column because it would slow calculation too much otherwise.
    // NOTICE: I add W when in last block because it is actually result from W cells to the left and W cells up.
    k = min(k, bl->score
    + max(targetLength - c - 1, queryLength - ((1 + lastBlock) * WORD_SIZE - 1) - 1)
      + (lastBlock == maxNumBlocks - 1 ? W : 0));
      
      //---------- Adjust number of blocks according to Ukkonen ----------//
      //--- Adjust last block ---//
      // If block is not beneath band, calculate next block. Only next because others are certainly beneath band.
      if (lastBlock + 1 < maxNumBlocks
      && !(//score[lastBlock] >= k + WORD_SIZE ||  // NOTICE: this condition could be satisfied if above block also!
          ((lastBlock + 1) * WORD_SIZE - 1
      > k - bl->score + 2 * WORD_SIZE - 2 - targetLength + c + queryLength))) {
        lastBlock++; bl++;
        bl->P = static_cast<Word>(-1); // All 1s
        bl->M = static_cast<Word>(0);
        int newHout = calculateBlock(bl->P, bl->M, Peq_c[lastBlock], hout, bl->P, bl->M);
        bl->score = (bl - 1)->score - hout + WORD_SIZE + newHout;
        hout = newHout;
      }
      
      // While block is out of band, move one block up.
      // NOTE: Condition used here is more loose than the one from the article, since I simplified the max() part of it.
      // I could consider adding that max part, for optimal performance.
      while (lastBlock >= firstBlock
      && (bl->score >= k + WORD_SIZE
      || ((lastBlock + 1) * WORD_SIZE - 1 >
      // TODO: Does not work if do not put +1! Why???
      k - bl->score + 2 * WORD_SIZE - 2 - targetLength + c + queryLength + 1))) {
        lastBlock--; bl--;
      }
      //-------------------------//
      
      //--- Adjust first block ---//
      // While outside of band, advance block
      while (firstBlock <= lastBlock
      && (blocks[firstBlock].score >= k + WORD_SIZE
      || ((firstBlock + 1) * WORD_SIZE - 1 <
        blocks[firstBlock].score - k - targetLength + queryLength + c))) {
        firstBlock++;
      }
      //--------------------------/
      
      
      // TODO: consider if this part is useful, it does not seem to help much
      if (c % STRONG_REDUCE_NUM == 0) { // Every some columns do more expensive but more efficient reduction
        while (lastBlock >= firstBlock) {
          // If all cells outside of band, remove block
          vector<int> scores = getBlockCellValues(*bl);
          int numCells = lastBlock == maxNumBlocks - 1 ? WORD_SIZE - W : WORD_SIZE;
          int r = lastBlock * WORD_SIZE + numCells - 1;
          bool reduce = true;
          for (int i = WORD_SIZE - numCells; i < WORD_SIZE; i++) {
            // TODO: Does not work if do not put +1! Why???
            if (scores[i] <= k && r <= k - scores[i] - targetLength + c + queryLength + 1) {
              reduce = false;
              break;
            }
            r--;
          }
          if (!reduce) break;
          lastBlock--; bl--;
        }
        
        while (firstBlock <= lastBlock) {
          // If all cells outside of band, remove block
          vector<int> scores = getBlockCellValues(blocks[firstBlock]);
          int numCells = firstBlock == maxNumBlocks - 1 ? WORD_SIZE - W : WORD_SIZE;
          int r = firstBlock * WORD_SIZE + numCells - 1;
          bool reduce = true;
          for (int i = WORD_SIZE - numCells; i < WORD_SIZE; i++) {
            if (scores[i] <= k && r >= scores[i] - k - targetLength + c + queryLength) {
              reduce = false;
              break;
            }
            r--;
          }
          if (!reduce) break;
          firstBlock++;
        }
      }
      
      
      // If band stops to exist finish
      if (lastBlock < firstBlock) {
        *bestScore_ = *position_ = -1;
        delete[] blocks;
        return EDLIB_STATUS_OK;
      }
      //------------------------------------------------------------------//
      
      
      //---- Save column so it can be used for reconstruction ----//
      if (findAlignment && c < targetLength) {
        bl = blocks + firstBlock;
        for (int b = firstBlock; b <= lastBlock; b++) {
          (*alignData)->Ps[maxNumBlocks * c + b] = bl->P;
          (*alignData)->Ms[maxNumBlocks * c + b] = bl->M;
          (*alignData)->scores[maxNumBlocks * c + b] = bl->score;
          bl++;
        }
        (*alignData)->firstBlocks[c] = firstBlock;
        (*alignData)->lastBlocks[c] = lastBlock;
      }
      //----------------------------------------------------------//
      //---- If this is stop column, save it and finish ----//
      if (c == targetStopPosition) {
        for (int b = firstBlock; b <= lastBlock; b++) {
          (*alignData)->Ps[b] = (blocks + b)->P;
          (*alignData)->Ms[b] = (blocks + b)->M;
          (*alignData)->scores[b] = (blocks + b)->score;
        }
        (*alignData)->firstBlocks[0] = firstBlock;
        (*alignData)->lastBlocks[0] = lastBlock;
        *bestScore_ = -1;
        *position_ = targetStopPosition;
        delete[] blocks;
        return EDLIB_STATUS_OK;
      }
      //----------------------------------------------------//
      
      targetChar++;
  }
  
  if (lastBlock == maxNumBlocks - 1) { // If last block of last column was calculated
    // Obtain best score from block -> it is complicated because query is padded with W cells
    int bestScore = getBlockCellValues(blocks[lastBlock])[W];
    if (bestScore <= k) {
      *bestScore_ = bestScore;
      *position_ = targetLength - 1;
      delete[] blocks;
      return EDLIB_STATUS_OK;
    }
  }
  
  *bestScore_ = *position_ = -1;
  delete[] blocks;
  return EDLIB_STATUS_OK;
}


/**
 * Finds one possible alignment that gives optimal score by moving back through the dynamic programming matrix,
 * that is stored in alignData. Consumes large amount of memory: O(queryLength * targetLength).
 * @param [in] queryLength  Normal length, without W.
 * @param [in] targetLength  Normal length, without W.
 * @param [in] bestScore  Best score.
 * @param [in] alignData  Data obtained during finding best score that is useful for finding alignment.
 * @param [out] alignment  Alignment.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignmentTraceback(const int queryLength, const int targetLength,
  const int bestScore, const AlignmentData* const alignData,
  unsigned char** const alignment, int* const alignmentLength) {
  const int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  const int W = maxNumBlocks * WORD_SIZE - queryLength;
  
  *alignment = static_cast<unsigned char*>(malloc((queryLength + targetLength - 1) * sizeof(unsigned char)));
  *alignmentLength = 0;
  int c = targetLength - 1; // index of column
  int b = maxNumBlocks - 1; // index of block in column
  int currScore = bestScore; // Score of current cell
  int lScore  = -1; // Score of left cell
  int uScore  = -1; // Score of upper cell
  int ulScore = -1; // Score of upper left cell
  Word currP = alignData->Ps[c * maxNumBlocks + b]; // P of current block
  Word currM = alignData->Ms[c * maxNumBlocks + b]; // M of current block
  // True if block to left exists and is in band
  bool thereIsLeftBlock = c > 0 && b >= alignData->firstBlocks[c-1] && b <= alignData->lastBlocks[c-1];
  // We set initial values of lP and lM to 0 only to avoid compiler warnings, they should not affect the
  // calculation as both lP and lM should be initialized at some moment later (but compiler can not
  // detect it since this initialization is guaranteed by "business" logic).
  Word lP = 0, lM = 0;
  if (thereIsLeftBlock) {
    lP = alignData->Ps[(c - 1) * maxNumBlocks + b]; // P of block to the left
    lM = alignData->Ms[(c - 1) * maxNumBlocks + b]; // M of block to the left
  }
  currP <<= W;
  currM <<= W;
  int blockPos = WORD_SIZE - W - 1; // 0 based index of current cell in blockPos
  
  // TODO(martin): refactor this whole piece of code. There are too many if-else statements,
  // it is too easy for a bug to hide and to hard to effectively cover all the edge-cases.
  // We need better separation of logic and responsibilities.
  while (true) {
    if (c == 0) {
      thereIsLeftBlock = true;
      lScore = b * WORD_SIZE + blockPos + 1;
      ulScore = lScore - 1;
    }
    
    // TODO: improvement: calculate only those cells that are needed,
    //       for example if I calculate upper cell and can move up,
    //       there is no need to calculate left and upper left cell
    //---------- Calculate scores ---------//
    if (lScore == -1 && thereIsLeftBlock) {
      lScore = alignData->scores[(c - 1) * maxNumBlocks + b]; // score of block to the left
      for (int i = 0; i < WORD_SIZE - blockPos - 1; i++) {
        if (lP & HIGH_BIT_MASK) lScore--;
        if (lM & HIGH_BIT_MASK) lScore++;
        lP <<= 1;
        lM <<= 1;
      }
    }
    if (ulScore == -1) {
      if (lScore != -1) {
        ulScore = lScore;
        if (lP & HIGH_BIT_MASK) ulScore--;
        if (lM & HIGH_BIT_MASK) ulScore++;
      }
      else if (c > 0 && b-1 >= alignData->firstBlocks[c-1] && b-1 <= alignData->lastBlocks[c-1]) {
        // This is the case when upper left cell is last cell in block,
        // and block to left is not in band so lScore is -1.
        ulScore = alignData->scores[(c - 1) * maxNumBlocks + b - 1];
      }
    }
    if (uScore == -1) {
      uScore = currScore;
      if (currP & HIGH_BIT_MASK) uScore--;
      if (currM & HIGH_BIT_MASK) uScore++;
      currP <<= 1;
      currM <<= 1;
    }
    //-------------------------------------//
    
    // TODO: should I check if there is upper block?
    
    //-------------- Move --------------//
    // Move up - insertion to target - deletion from query
    if (uScore != -1 && uScore + 1 == currScore) {
      currScore = uScore;
      lScore = ulScore;
      uScore = ulScore = -1;
      if (blockPos == 0) { // If entering new (upper) block
        if (b == 0) { // If there are no cells above (only boundary cells)
          (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT; // Move up
          for (int i = 0; i < c + 1; i++) // Move left until end
            (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
          break;
        } else {
          blockPos = WORD_SIZE - 1;
          b--;
          currP = alignData->Ps[c * maxNumBlocks + b];
          currM = alignData->Ms[c * maxNumBlocks + b];
          if (c > 0 && b >= alignData->firstBlocks[c-1] && b <= alignData->lastBlocks[c-1]) {
            thereIsLeftBlock = true;
            lP = alignData->Ps[(c - 1) * maxNumBlocks + b]; // TODO: improve this, too many operations
            lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
          } else {
            thereIsLeftBlock = false;
            // TODO(martin): There may not be left block, but there can be left boundary - do we
            // handle this correctly then? Are l and ul score set correctly? I should check that / refactor this.
          }
        }
      } else {
        blockPos--;
        lP <<= 1;
        lM <<= 1;
      }
      // Mark move
      (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
    }
    // Move left - deletion from target - insertion to query
    else if (lScore != -1 && lScore + 1 == currScore) {
      currScore = lScore;
      uScore = ulScore;
      lScore = ulScore = -1;
      c--;
      if (c == -1) { // If there are no cells to the left (only boundary cells)
        (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE; // Move left
        int numUp = b * WORD_SIZE + blockPos + 1;
        for (int i = 0; i < numUp; i++) // Move up until end
          (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
        break;
      }
      currP = lP;
      currM = lM;
      if (c > 0 && b >= alignData->firstBlocks[c-1] && b <= alignData->lastBlocks[c-1]) {
        thereIsLeftBlock = true;
        lP = alignData->Ps[(c - 1) * maxNumBlocks + b];
        lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
      } else {
        if (c == 0) { // If there are no cells to the left (only boundary cells)
          thereIsLeftBlock = true;
          lScore = b * WORD_SIZE + blockPos + 1;
          ulScore = lScore - 1;
        } else {
          thereIsLeftBlock = false;
        }
      }
      // Mark move
      (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
    }
    // Move up left - (mis)match
    else if (ulScore != -1) {
      unsigned char moveCode = ulScore == currScore ? EDLIB_EDOP_MATCH : EDLIB_EDOP_MISMATCH;
      currScore = ulScore;
      uScore = lScore = ulScore = -1;
      c--;
      if (c == -1) { // If there are no cells to the left (only boundary cells)
        (*alignment)[(*alignmentLength)++] = moveCode; // Move left
        int numUp = b * WORD_SIZE + blockPos;
        for (int i = 0; i < numUp; i++) // Move up until end
          (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
        break;
      }
      if (blockPos == 0) { // If entering upper left block
        if (b == 0) { // If there are no more cells above (only boundary cells)
          (*alignment)[(*alignmentLength)++] = moveCode; // Move up left
          for (int i = 0; i < c + 1; i++) // Move left until end
            (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
          break;
        }
        blockPos = WORD_SIZE - 1;
        b--;
        currP = alignData->Ps[c * maxNumBlocks + b];
        currM = alignData->Ms[c * maxNumBlocks + b];
      } else { // If entering left block
        blockPos--;
        currP = lP;
        currM = lM;
        currP <<= 1;
        currM <<= 1;
      }
      // Set new left block
      if (c > 0 && b >= alignData->firstBlocks[c-1] && b <= alignData->lastBlocks[c-1]) {
        thereIsLeftBlock = true;
        lP = alignData->Ps[(c - 1) * maxNumBlocks + b];
        lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
      } else {
        if (c == 0) { // If there are no cells to the left (only boundary cells)
          thereIsLeftBlock = true;
          lScore = b * WORD_SIZE + blockPos + 1;
          ulScore = lScore - 1;
        } else {
          thereIsLeftBlock = false;
        }
      }
      // Mark move
      (*alignment)[(*alignmentLength)++] = moveCode;
    } else {
      // Reached end - finished!
      break;
    }
    //----------------------------------//
  }
  
  *alignment = static_cast<unsigned char*>(realloc(*alignment, (*alignmentLength) * sizeof(unsigned char)));
  reverse(*alignment, *alignment + (*alignmentLength));
  return EDLIB_STATUS_OK;
}


/**
 * Finds one possible alignment that gives optimal score (bestScore).
 * It will split problem into smaller problems using Hirschberg's algorithm and when they are small enough,
 * it will solve them using traceback algorithm.
 * @param [in] query
 * @param [in] rQuery  Reversed query.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] rTarget  Reversed target.
 * @param [in] targetLength
 * @param [in] equalityDefinition
 * @param [in] alphabetLength
 * @param [in] bestScore  Best(optimal) score.
 * @param [out] alignment  Sequence of edit operations that make target equal to query.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignment(
    const unsigned char* const query, const unsigned char* const rQuery, const int queryLength,
    const unsigned char* const target, const unsigned char* const rTarget, const int targetLength,
    const EqualityDefinition& equalityDefinition, const int alphabetLength, const int bestScore,
    unsigned char** const alignment, int* const alignmentLength) {
  
  // Handle special case when one of sequences has length of 0.
  if (queryLength == 0 || targetLength == 0) {
    *alignmentLength = targetLength + queryLength;
    *alignment = static_cast<unsigned char*>(malloc((*alignmentLength) * sizeof(unsigned char)));
    for (int i = 0; i < *alignmentLength; i++) {
      (*alignment)[i] = queryLength == 0 ? EDLIB_EDOP_DELETE : EDLIB_EDOP_INSERT;
    }
    return EDLIB_STATUS_OK;
  }
  
  const int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  const int W = maxNumBlocks * WORD_SIZE - queryLength;
  int statusCode;
  
  // TODO: think about reducing number of memory allocations in alignment functions, probably
  // by sharing some memory that is allocated only once. That refers to: Peq, columns in Hirschberg,
  // and it could also be done for alignments - we could have one big array for alignment that would be
  // sparsely populated by each of steps in recursion, and at the end we would just consolidate those results.
  
  // If estimated memory consumption for traceback algorithm is smaller than 1MB use it,
  // otherwise use Hirschberg's algorithm. By running few tests I choose boundary of 1MB as optimal.
  long long alignmentDataSize = (2ll * sizeof(Word) + sizeof(int)) * maxNumBlocks * targetLength
  + 2ll * sizeof(int) * targetLength;
  if (alignmentDataSize < 1024 * 1024) {
    int score_, endLocation_;  // Used only to call function.
    AlignmentData* alignData = NULL;
    Word* Peq = buildPeq(alphabetLength, query, queryLength, equalityDefinition);
    myersCalcEditDistanceNW(Peq, W, maxNumBlocks,
      queryLength,
      target, targetLength,
      bestScore,
      &score_, &endLocation_, true, &alignData, -1);
    //assert(score_ == bestScore);
    //assert(endLocation_ == targetLength - 1);
    
    statusCode = obtainAlignmentTraceback(queryLength, targetLength,
      bestScore, alignData, alignment, alignmentLength);
    delete alignData;
    delete[] Peq;
  } else {
    statusCode = obtainAlignmentHirschberg(query, rQuery, queryLength,
      target, rTarget, targetLength,
      equalityDefinition, alphabetLength, bestScore,
      alignment, alignmentLength);
  }
  return statusCode;
}


/**
 * Finds one possible alignment that gives optimal score (bestScore).
 * Uses Hirschberg's algorithm to split problem into two sub-problems, solve them and combine them together.
 * @param [in] query
 * @param [in] rQuery  Reversed query.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] rTarget  Reversed target.
 * @param [in] targetLength
 * @param [in] alphabetLength
 * @param [in] bestScore  Best(optimal) score.
 * @param [out] alignment  Sequence of edit operations that make target equal to query.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignmentHirschberg(
    const unsigned char* const query, const unsigned char* const rQuery, const int queryLength,
    const unsigned char* const target, const unsigned char* const rTarget, const int targetLength,
    const EqualityDefinition& equalityDefinition, const int alphabetLength, const int bestScore,
    unsigned char** const alignment, int* const alignmentLength) {
  
  const int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
  const int W = maxNumBlocks * WORD_SIZE - queryLength;
  
  Word* Peq = buildPeq(alphabetLength, query, queryLength, equalityDefinition);
  Word* rPeq = buildPeq(alphabetLength, rQuery, queryLength, equalityDefinition);
  
  // Used only to call functions.
  int score_, endLocation_;
  
  // Divide dynamic matrix into two halfs, left and right.
  const int leftHalfWidth = targetLength / 2;
  const int rightHalfWidth = targetLength - leftHalfWidth;
  
  // Calculate left half.
  AlignmentData* alignDataLeftHalf = NULL;
  int leftHalfCalcStatus = myersCalcEditDistanceNW(
    Peq, W, maxNumBlocks, queryLength, target, targetLength, bestScore,
    &score_, &endLocation_, false, &alignDataLeftHalf, leftHalfWidth - 1);
  
  // Calculate right half.
  AlignmentData* alignDataRightHalf = NULL;
  int rightHalfCalcStatus = myersCalcEditDistanceNW(
    rPeq, W, maxNumBlocks, queryLength, rTarget, targetLength, bestScore,
    &score_, &endLocation_, false, &alignDataRightHalf, rightHalfWidth - 1);
  
  delete[] Peq;
  delete[] rPeq;
  
  if (leftHalfCalcStatus == EDLIB_STATUS_ERROR || rightHalfCalcStatus == EDLIB_STATUS_ERROR) {
    if (alignDataLeftHalf) delete alignDataLeftHalf;
    if (alignDataRightHalf) delete alignDataRightHalf;
    return EDLIB_STATUS_ERROR;
  }
  
  // Unwrap the left half.
  int firstBlockIdxLeft = alignDataLeftHalf->firstBlocks[0];
  int lastBlockIdxLeft = alignDataLeftHalf->lastBlocks[0];
  // TODO: avoid this allocation by using some shared array?
  // scoresLeft contains scores from left column, starting with scoresLeftStartIdx row (query index)
  // and ending with scoresLeftEndIdx row (0-indexed).
  int scoresLeftLength = (lastBlockIdxLeft - firstBlockIdxLeft + 1) * WORD_SIZE;
  int* scoresLeft = new int[scoresLeftLength];
  for (int blockIdx = firstBlockIdxLeft; blockIdx <= lastBlockIdxLeft; blockIdx++) {
    Block block(alignDataLeftHalf->Ps[blockIdx], alignDataLeftHalf->Ms[blockIdx],
      alignDataLeftHalf->scores[blockIdx]);
    readBlock(block, scoresLeft + (blockIdx - firstBlockIdxLeft) * WORD_SIZE);
  }
  int scoresLeftStartIdx = firstBlockIdxLeft * WORD_SIZE;
  // If last block contains padding, shorten the length of scores for the length of padding.
  if (lastBlockIdxLeft == maxNumBlocks - 1) {
    scoresLeftLength -= W;
  }
  
  // Unwrap the right half (I also reverse it while unwraping).
  int firstBlockIdxRight = alignDataRightHalf->firstBlocks[0];
  int lastBlockIdxRight = alignDataRightHalf->lastBlocks[0];
  int scoresRightLength = (lastBlockIdxRight - firstBlockIdxRight + 1) * WORD_SIZE;
  int* scoresRight = new int[scoresRightLength];
  int* scoresRightOriginalStart = scoresRight;
  for (int blockIdx = firstBlockIdxRight; blockIdx <= lastBlockIdxRight; blockIdx++) {
    Block block(alignDataRightHalf->Ps[blockIdx], alignDataRightHalf->Ms[blockIdx],
      alignDataRightHalf->scores[blockIdx]);
    readBlockReverse(block, scoresRight + (lastBlockIdxRight - blockIdx) * WORD_SIZE);
  }
  int scoresRightStartIdx = queryLength - (lastBlockIdxRight + 1) * WORD_SIZE;
  // If there is padding at the beginning of scoresRight (that can happen because of reversing that we do),
  // move pointer forward to remove the padding (that is why we remember originalStart).
  if (scoresRightStartIdx < 0) {
    //assert(scoresRightStartIdx == -1 * W);
    scoresRight += W;
    scoresRightStartIdx += W;
    scoresRightLength -= W;
  }
  
  delete alignDataLeftHalf;
  delete alignDataRightHalf;
  
  //--------------------- Find the best move ----------------//
  // Find the query/row index of cell in left column which together with its lower right neighbour
  // from right column gives the best score (when summed). We also have to consider boundary cells
  // (those cells at -1 indexes).
  //  x|
  //  -+-
  //   |x
  int queryIdxLeftStart = max(scoresLeftStartIdx, scoresRightStartIdx - 1);
  int queryIdxLeftEnd = min(scoresLeftStartIdx + scoresLeftLength - 1,
    scoresRightStartIdx + scoresRightLength - 2);
  int leftScore = -1, rightScore = -1;
  int queryIdxLeftAlignment = -1;  // Query/row index of cell in left column where alignment is passing through.
  bool queryIdxLeftAlignmentFound = false;
  for (int queryIdx = queryIdxLeftStart; queryIdx <= queryIdxLeftEnd; queryIdx++) {
    leftScore = scoresLeft[queryIdx - scoresLeftStartIdx];
    rightScore = scoresRight[queryIdx + 1 - scoresRightStartIdx];
    if (leftScore + rightScore == bestScore) {
      queryIdxLeftAlignment = queryIdx;
      queryIdxLeftAlignmentFound = true;
      break;
    }
  }
  // Check boundary cells.
  if (!queryIdxLeftAlignmentFound && scoresLeftStartIdx == 0 && scoresRightStartIdx == 0) {
    leftScore = leftHalfWidth;
    rightScore = scoresRight[0];
    if (leftScore + rightScore == bestScore) {
      queryIdxLeftAlignment = -1;
      queryIdxLeftAlignmentFound = true;
    }
  }
  if (!queryIdxLeftAlignmentFound && scoresLeftStartIdx + scoresLeftLength == queryLength
  && scoresRightStartIdx + scoresRightLength == queryLength) {
    leftScore = scoresLeft[scoresLeftLength - 1];
    rightScore = rightHalfWidth;
    if (leftScore + rightScore == bestScore) {
      queryIdxLeftAlignment = queryLength - 1;
      queryIdxLeftAlignmentFound = true;
    }
  }
  
  delete[] scoresLeft;
  delete[] scoresRightOriginalStart;
  
  if (queryIdxLeftAlignmentFound == false) {
    // If there was no move that is part of optimal alignment, then there is no such alignment
    // or given bestScore is not correct!
    return EDLIB_STATUS_ERROR;
  }
  //----------------------------------------------------------//
  
  // Calculate alignments for upper half of left half (upper left - ul)
  // and lower half of right half (lower right - lr).
  const int ulHeight = queryIdxLeftAlignment + 1;
  const int lrHeight = queryLength - ulHeight;
  const int ulWidth = leftHalfWidth;
  const int lrWidth = rightHalfWidth;
  unsigned char* ulAlignment = NULL; int ulAlignmentLength;
  int ulStatusCode = obtainAlignment(query, rQuery + lrHeight, ulHeight,
    target, rTarget + lrWidth, ulWidth,
    equalityDefinition, alphabetLength, leftScore,
    &ulAlignment, &ulAlignmentLength);
  unsigned char* lrAlignment = NULL; int lrAlignmentLength;
  int lrStatusCode = obtainAlignment(query + ulHeight, rQuery, lrHeight,
    target + ulWidth, rTarget, lrWidth,
    equalityDefinition, alphabetLength, rightScore,
    &lrAlignment, &lrAlignmentLength);
  if (ulStatusCode == EDLIB_STATUS_ERROR || lrStatusCode == EDLIB_STATUS_ERROR) {
    if (ulAlignment) free(ulAlignment);
    if (lrAlignment) free(lrAlignment);
    return EDLIB_STATUS_ERROR;
  }
  
  // Build alignment by concatenating upper left alignment with lower right alignment.
  *alignmentLength = ulAlignmentLength + lrAlignmentLength;
  *alignment = static_cast<unsigned char*>(malloc((*alignmentLength) * sizeof(unsigned char)));
  memcpy(*alignment, ulAlignment, ulAlignmentLength);
  memcpy(*alignment + ulAlignmentLength, lrAlignment, lrAlignmentLength);
  
  free(ulAlignment);
  free(lrAlignment);
  return EDLIB_STATUS_OK;
}


/**
 * Takes char query and char target, recognizes alphabet and transforms them into unsigned char sequences
 * where elements in sequences are not any more letters of alphabet, but their index in alphabet.
 * Most of internal edlib functions expect such transformed sequences.
 * This function will allocate queryTransformed and targetTransformed, so make sure to free them when done.
 * Example:
 *   Original sequences: "ACT" and "CGT".
 *   Alphabet would be recognized as "ACTG". Alphabet length = 4.
 *   Transformed sequences: [0, 1, 2] and [1, 3, 2].
 * @param [in] queryOriginal
 * @param [in] queryLength
 * @param [in] targetOriginal
 * @param [in] targetLength
 * @param [out] queryTransformed  It will contain values in range [0, alphabet length - 1].
 * @param [out] targetTransformed  It will contain values in range [0, alphabet length - 1].
 * @return  Alphabet as a string of unique characters, where index of each character is its value in transformed
 *          sequences.
 */
static string transformSequences(const char* const queryOriginal, const int queryLength,
  const char* const targetOriginal, const int targetLength,
  unsigned char** const queryTransformed,
  unsigned char** const targetTransformed) {
  // Alphabet is constructed from letters that are present in sequences.
  // Each letter is assigned an ordinal number, starting from 0 up to alphabetLength - 1,
  // and new query and target are created in which letters are replaced with their ordinal numbers.
  // This query and target are used in all the calculations later.
  *queryTransformed = static_cast<unsigned char *>(malloc(sizeof(unsigned char) * queryLength));
  *targetTransformed = static_cast<unsigned char *>(malloc(sizeof(unsigned char) * targetLength));
  
  string alphabet = "";
  
  // Alphabet information, it is constructed on fly while transforming sequences.
  // letterIdx[c] is index of letter c in alphabet.
  unsigned char letterIdx[MAX_UCHAR + 1];
  bool inAlphabet[MAX_UCHAR + 1]; // inAlphabet[c] is true if c is in alphabet
  for (int i = 0; i < MAX_UCHAR + 1; i++) inAlphabet[i] = false;
  
  for (int i = 0; i < queryLength; i++) {
    unsigned char c = static_cast<unsigned char>(queryOriginal[i]);
    if (!inAlphabet[c]) {
      inAlphabet[c] = true;
      letterIdx[c] = static_cast<unsigned char>(alphabet.size());
      alphabet += queryOriginal[i];
    }
    (*queryTransformed)[i] = letterIdx[c];
  }
  for (int i = 0; i < targetLength; i++) {
    unsigned char c = static_cast<unsigned char>(targetOriginal[i]);
    if (!inAlphabet[c]) {
      inAlphabet[c] = true;
      letterIdx[c] = static_cast<unsigned char>(alphabet.size());
      alphabet += targetOriginal[i];
    }
    (*targetTransformed)[i] = letterIdx[c];
  }
  
  return alphabet;
}


extern "C" EdlibAlignConfig edlibNewAlignConfig(int k, EdlibAlignMode mode, EdlibAlignTask task,
  const EdlibEqualityPair* additionalEqualities,
  int additionalEqualitiesLength) {
  EdlibAlignConfig config;
  config.k = k;
  config.mode = mode;
  config.task = task;
  config.additionalEqualities = additionalEqualities;
  config.additionalEqualitiesLength = additionalEqualitiesLength;
  return config;
}

extern "C" EdlibAlignConfig edlibDefaultAlignConfig(void) {
  return edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
}

extern "C" void edlibFreeAlignResult(EdlibAlignResult result) {
  if (result.endLocations) free(result.endLocations);
  if (result.startLocations) free(result.startLocations);
  if (result.alignment) free(result.alignment);
}

// edlibR str_length input
std::size_t strlen_utf8(const std::string& str) {
  std::size_t length = 0;
  for (char c : str) {
    if ((c & 0xC0) != 0x80) {
      ++length;
    }
  }
  return length;
}

std::string revcomp(const std::string& sequence) {
  std::unordered_map<char, char> complement {
    {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'},
    {'a', 't'}, {'t', 'a'}, {'c', 'g'}, {'g', 'c'},
    {'N', 'N'}, {'n', 'n'}
  };
  std::string reversed(sequence.rbegin(), sequence.rend());
  for (char& nucleotide : reversed) {
    nucleotide = complement[nucleotide];
  }
  return reversed;
}

int dna_hamming_distance(const std::string& a, const std::string& b) {
  int count = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] != b[i]) {
      ++count;
    }
  }
  return count;
}

std::vector<int64_t> barcodes_to_bits_cpp(const std::vector<std::string>& barcodes) {
  int n = barcodes.size();
  std::vector<int64_t> results(n);
  for(int i = 0; i < n; ++i) {
    const std::string& barcode = barcodes[i];
    int64_t result = 0;
    for (char c : barcode) {
      result <<= 2;
      switch (c) {
      case 'A': break;
      case 'C': result |= 1; break;
      case 'T': result |= 2; break;
      case 'G': result |= 3; break;
      }
    }
    results[i] = result;
  }
  return results;
}

std::vector<std::string> bits_to_barcodes_cpp(const std::vector<int64_t>& input, int barcode_length = 16) {
  int n = input.size();
  std::vector<std::string> results(n);
  for(int i = 0; i < n; ++i) {
    int64_t int64_code = input[i];
    std::string result;
    for (int j = 0; j < barcode_length; ++j) {
      switch (int64_code & 3) {
      case 0: result.insert(result.begin(), 'A'); break;
      case 1: result.insert(result.begin(), 'C'); break;
      case 2: result.insert(result.begin(), 'T'); break;
      case 3: result.insert(result.begin(), 'G'); break;
      }
      int64_code >>= 2;
    }
    results[i] = result;
  }
  return results;
}

void recursive_flip_bits_cpp(int64_t original_value, int64_t value, int depth, int max_iterations,
  int num_bits, const std::unordered_set<int64_t>& whitelist,
  int& min_hamming_distance, std::set<int64_t>& min_hamming_results) {
  if (depth >= max_iterations) return;
  for (int index = 0; index < num_bits; ++index) {
    int64_t mask = 1LL << (num_bits - 1 - index);
    int64_t result = value ^ mask;
    if (result != original_value) {
      if (whitelist.find(result) != whitelist.end()) {
        int hamming_dist = dna_hamming_distance(bits_to_barcodes_cpp({original_value})[0], bits_to_barcodes_cpp({result})[0]);
        if (hamming_dist < min_hamming_distance) {
          min_hamming_distance = hamming_dist;
          min_hamming_results.clear();
          min_hamming_results.insert(result);
        } else if (hamming_dist == min_hamming_distance) {
          min_hamming_results.insert(result);
        }
      }
      recursive_flip_bits_cpp(original_value, result, depth + 1, max_iterations, num_bits, whitelist, min_hamming_distance, min_hamming_results);
    }
  }
}

std::map<int64_t, std::vector<std::pair<int64_t, int>>> 
  mutate_and_check_list_cpp(const std::vector<int64_t>& barcodes, 
    const std::unordered_set<int64_t>& whitelist, 
    int max_iterations = 4, int max_mutations = 4, 
    int barcode_length = 16) {
    std::map<int64_t, std::vector<std::pair<int64_t, int>>> results;
    int num_bits = 2 * barcode_length;
    for (const auto& int64_barcode : barcodes) {
      int min_hamming_distance = max_mutations; // Initialize with max_mutations
      std::set<int64_t> min_hamming_results;
      recursive_flip_bits_cpp(int64_barcode, int64_barcode, 0, max_iterations, num_bits, 
        whitelist, min_hamming_distance, min_hamming_results);
      if (!min_hamming_results.empty()) {
        std::vector<std::pair<int64_t, int>> mutated_barcodes;
        for (const auto& mutated : min_hamming_results) {
          mutated_barcodes.emplace_back(mutated, min_hamming_distance);
        }
        results[int64_barcode] = mutated_barcodes;
      }
    }
    return results;
  }

void recursive_flip_bits_list(int64_t original_value, int64_t value, int depth, int max_mutations,
  int num_bits, std::set<int64_t>& results, std::unordered_set<int64_t>& whitelist) {
  if (depth >= max_mutations) return;
  for (int index = 0; index < num_bits; ++index) {
    int64_t mask = 1LL << (num_bits - 1 - index);
    int64_t result = value ^ mask;
    if (result != original_value) {
      if (whitelist.find(result) != whitelist.end()) {
        results.insert(result);
      }
      recursive_flip_bits_list(original_value, result, depth + 1, max_mutations, num_bits, results, whitelist);
    }
  }
}

std::unordered_map<std::string, std::vector<int>> sig_parse(const std::string& signature) {
  std::unordered_map<std::string, std::vector<int>> parsedSignature;
  std::stringstream ss(signature);
  std::string token;
  while (std::getline(ss, token, '|')) {
    size_t pos = token.find(":");
    std::string id = token.substr(0, pos);
    std::string coords = token.substr(pos + 1);
    std::stringstream coord_ss(coords);
    std::string coord;
    std::vector<int> coord_vec;
    while (std::getline(coord_ss, coord, ':')) {
      coord_vec.push_back(std::stoi(coord));
    }
    parsedSignature[id] = coord_vec;
  }
  return parsedSignature;
}

std::vector<int> find_mismatch_positions(const std::string& original, const std::string& mutated) {
  std::vector<int> mismatch_positions;
  for (size_t i = 0; i < original.size(); ++i) {
    if (original[i] != mutated[i]) {
      mismatch_positions.push_back(i + 1);  // Add 1 to convert from 0-based to 1-based index
    }
  }
  return mismatch_positions;
}

std::pair<std::vector<std::string>, double> best_matches_jaccard_distance_qgram(const std::string& original, const std::vector<std::string>& candidates, int q) {
  double min_jaccard_distance = 1.0;
  std::vector<std::string> best_matches;
  for (const auto& candidate : candidates) {
    std::unordered_set<std::string> setA, setB;
    for (int i = 0; i <= static_cast<int>(original.size()) - q; ++i) {
      setA.insert(original.substr(i, q));
    }
    for (int i = 0; i <= static_cast<int>(candidate.size()) - q; ++i) {
      setB.insert(candidate.substr(i, q));
    }
    std::unordered_set<std::string> intersection;
    for (const auto& item : setA) {
      if (setB.find(item) != setB.end()) {
        intersection.insert(item);
      }
    }
    double jaccard_index = static_cast<double>(intersection.size()) / (setA.size() + setB.size() - intersection.size());
    double jaccard_distance = 1 - jaccard_index;
    if (jaccard_distance < min_jaccard_distance) {
      min_jaccard_distance = jaccard_distance;
      best_matches.clear();
      best_matches.push_back(candidate);
    } else if (jaccard_distance == min_jaccard_distance) {
      best_matches.push_back(candidate);
    }
  }
  return std::make_pair(best_matches, min_jaccard_distance);
}

std::vector<std::string> split(const std::string& s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

std::unordered_map<std::string, double> calculateLowerBounds(DataFrame null_distance) {
  std::unordered_map<std::string, double> lowerBounds;
  CharacterVector query_id = null_distance["query_id"];
  NumericVector null_dist = null_distance["null_distance"];
  NumericVector sd_null = null_distance["sd_null"];
  for (int i = 0; i < query_id.size(); ++i) {
    lowerBounds[as<std::string>(query_id[i])] = std::floor(null_dist[i] - std::ceil(sd_null[i]));
  } 
  return lowerBounds;
}

std::unordered_map<std::string, 
  std::unordered_map<std::string,
    std::tuple<int, int, std::string>>> convertReadLayout(Rcpp::DataFrame read_layout) {
      std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<int, int, std::string>>> layoutMap;
      Rcpp::StringVector id = read_layout["id"];
      Rcpp::IntegerVector expected_length = read_layout["expected_length"];
      Rcpp::StringVector type = read_layout["type"];
      Rcpp::IntegerVector order = read_layout["order"];
      Rcpp::StringVector direction = read_layout["direction"];
      for (int i = 0; i < id.size(); ++i) {
        std::string id_str = Rcpp::as<std::string>(id[i]);
        int exp_len = expected_length[i];
        std::string type_str = Rcpp::as<std::string>(type[i]);
        int order_int = order[i];
        std::string direction_str = Rcpp::as<std::string>(direction[i]);
        layoutMap[id_str][direction_str] = std::make_tuple(exp_len, order_int, type_str);
      }
      return layoutMap;
    }

std::unordered_map<std::string, int> generateOrderMap(Rcpp::DataFrame read_layout) {
  std::unordered_map<std::string, int> orderMap;
  Rcpp::CharacterVector ids = read_layout["id"];
  Rcpp::IntegerVector orders = read_layout["order"];
  for (int i = 0; i < ids.size(); ++i) {
    std::string id = Rcpp::as<std::string>(ids[i]);
    int order = orders[i];
    orderMap[id] = order;
  }
  return orderMap;
}

std::unordered_map<std::string, int> countIDs(const std::unordered_map<std::string, std::vector<int>>& parsedSignature) {
  std::unordered_map<std::string, int> idCount;
  for (const auto& item : parsedSignature) {
    std::string id = item.first;
    const std::vector<int>& coords = item.second;
    idCount[id] += 1;
  } 
  return idCount;
} 

std::string error_categorization(const std::string& barcoded_sigstring, 
  const std::unordered_map<std::string, double>& lowerBounds,
  const std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<int, int, std::string>>>& layoutMap, 
  const std::unordered_map<std::string, int>& orderMap,
  bool verbose = false) {
  std::unordered_map<std::string, std::vector<int>> parsedSignature = sig_parse(barcoded_sigstring);
  if(verbose){
    List out;
    for (const auto& line : parsedSignature) {
      out[line.first] = line.second;
    }
    Rcpp::Rcout << "ParsedSignature output: ";
    Rf_PrintValue(out);
  }
  std::unordered_map<std::string, int> idCount = countIDs(parsedSignature);
  if(verbose){
    List out_id;
    for (const auto& line : idCount) {
      out_id[line.first] = line.second;
    }
    Rcpp::Rcout << "idCount output: ";
    Rf_PrintValue(out_id);
    Rcpp::Rcout << "Lower bounds at forw_primer: " << lowerBounds.at("forw_primer") << std::endl;
    Rcpp::Rcout << "Lower bounds at rc_forw_primer: " << lowerBounds.at("rc_forw_primer") << std::endl;
  }
  std::string category = "undecided";
  if (
      (
          // Existing conditions
          (parsedSignature.find("rc_forw_primer") == parsedSignature.end() || parsedSignature.at("rc_forw_primer")[0] >= lowerBounds.at("rc_forw_primer")) &&
            parsedSignature.at("forw_primer")[0] <= lowerBounds.at("forw_primer") &&
            parsedSignature.at("forw_primer")[0] < parsedSignature.at("rc_forw_primer")[0]
      ) ||
        ( // New condition: both primers are below the threshold
            parsedSignature.at("forw_primer")[0] <= lowerBounds.at("forw_primer") &&
              parsedSignature.at("rc_forw_primer")[0] <= lowerBounds.at("rc_forw_primer") &&
              (
                  parsedSignature.find("rev_primer") == parsedSignature.end() ||
                    parsedSignature.at("rev_primer")[0] <= lowerBounds.at("rev_primer")
              ) &&
                parsedSignature.at("forw_primer")[0] < parsedSignature.at("rev_primer")[0]
        )
  ){
    category = "forward";
    if(verbose){
      Rcpp::Rcout << "Mapped to the forward category, now attempting to find read coordinates."<< std::endl;
    }
    std::string read_key = "read";
    int read_order = orderMap.at(read_key);
    if(verbose){
      Rcpp::Rcout << "Now mapping the following sigstring: " << barcoded_sigstring << std::endl;
      Rcpp::Rcout << "Printing read_key and read order: " << read_order << " is read order & " << read_key << " is read key!" <<std::endl;
    }
    int new_start = -1;
    int new_end = -1;
    int search_distance = 1;
    std::vector<int> new_read_coords = {0, 0, 0};
    bool found_start = false, found_end = false;
    while (!found_start || !found_end) {
      if(verbose) {
        Rcpp::Rcout << "Current search distance: " << search_distance << std::endl;
      }
      for (const auto& item : parsedSignature) {
        std::string id = item.first;
        if (layoutMap.find(id) != layoutMap.end() && layoutMap.at(id).begin()->first == "forward") {
          int distance = abs(orderMap.at(id) - read_order);
          if (distance == search_distance) {
            if(verbose) {
              Rcpp::Rcout << "Checking adjacent ID: " << id << " at distance: " << distance << std::endl;
            }
            if (orderMap.at(id) < read_order) {
              if (verbose) {
                Rcpp::Rcout << "Checking id: " << id << std::endl;
              }
              if ((lowerBounds.find(id) != lowerBounds.end() && item.second[0] <= lowerBounds.at(id)) || id == "poly_a" || id == "poly_t" || id == "barcode" || id == "umi") {
                if (!found_start) {
                  new_start = item.second[2] + 1;
                }
                found_start = true; // Make sure to update the flag
                if(verbose) {
                  Rcpp::Rcout << "New start found: " << new_start << std::endl;
                }
                continue;
              }
            } else {
              if ((lowerBounds.find(id) != lowerBounds.end() && item.second[0] <= lowerBounds.at(id)) || id == "poly_a" || id == "poly_t" || id == "barcode" || id == "umi") {
                if(!found_end){
                  new_end = item.second[1] - 1;
                }
                found_end = true;
                if(verbose) {
                  Rcpp::Rcout << "New end found: " << new_end << std::endl;
                }
              }
              continue;
            }
          }
        }
      }
      if (new_end == -1 && new_start != -1) {
        if(verbose){
          Rcpp::Rcout << "New end not found, setting to sentinel value 0" << std::endl;
        }
        new_end = 0;  // Note the corrected assignment operator
        break;
      }
      if (found_start && found_end) {
        if(verbose) {
          Rcpp::Rcout << "Both new start and end found, exiting loop." << std::endl;
        }
        break;
      }
      if (!found_start || !found_end) {
        search_distance++;
        if(verbose) {
          Rcpp::Rcout << "New position not found, increasing search distance to: " << search_distance << std::endl;
        }
      }
    }
    if (new_start != -1 && new_end != -1) {
      new_read_coords = {0, new_start, new_end};
      if(verbose){
        Rcpp::Rcout << "putative new coordinates are " << new_read_coords[0] << new_read_coords[1] << new_read_coords[2] << std::endl;
      }
      parsedSignature[read_key] = new_read_coords;
    }
    std::stringstream newSignature;
    for (const auto& item : parsedSignature) {
      newSignature << "|" << item.first << ":" << item.second[0] << ":" << item.second[1] << ":" << item.second[2];
    }
    newSignature << "<" << category << ">";
    if(verbose){
      Rcpp::Rcout << "putative new signature is " << newSignature.str() <<  "!" << std::endl;
    }
    return newSignature.str().substr(1);
  }
  // Check for "reverses"
  if (
      (
          // Existing conditions
          (parsedSignature.find("forw_primer") == parsedSignature.end() || parsedSignature.at("forw_primer")[0] >= lowerBounds.at("forw_primer")) &&
            parsedSignature.at("rc_forw_primer")[0] <= lowerBounds.at("rc_forw_primer") &&
            parsedSignature.at("rc_forw_primer")[0] < parsedSignature.at("forw_primer")[0]
      ) ||
        (
            // New condition: both primers are below the threshold
            parsedSignature.at("forw_primer")[0] <= lowerBounds.at("forw_primer") &&
              parsedSignature.at("rc_forw_primer")[0] <= lowerBounds.at("rc_forw_primer") &&
              (
                  parsedSignature.find("rc_rev_primer") == parsedSignature.end() ||
                    parsedSignature.at("rc_rev_primer")[0] <= lowerBounds.at("rc_rev_primer")
              ) &&
                parsedSignature.at("rc_forw_primer")[0] < parsedSignature.at("rc_rev_primer")[0]
        )
  ){
    category = "reverse";
    if(verbose){
      Rcpp::Rcout << "Mapped to the reverse category, now attempting to find read coordinates."<< std::endl;
    }
    std::vector<int> new_read_coords = {0, 0, 0};
    std::string read_key = "rc_read";
    int read_order = orderMap.at(read_key);
    if(verbose){
      Rcpp::Rcout << "Now mapping the following sigstring: " << barcoded_sigstring << std::endl;
      Rcpp::Rcout << "Printing read_key and read order: " << read_order << " is read order & " << read_key << " is read key!" <<std::endl;
    }
    int new_start = -1;
    int new_end = -1;
    int search_distance = 1;
    bool found_start = false, found_end = false;
    while (!found_start || !found_end) {
      if(verbose) {
        Rcpp::Rcout << "Current search distance: " << search_distance << std::endl;
      }
      for (const auto& item : parsedSignature) {
        std::string id = item.first;
        if (layoutMap.find(id) != layoutMap.end() && layoutMap.at(id).begin()->first == "reverse") {
          int distance = abs(orderMap.at(id) - read_order);
          if (distance == search_distance) {
            if(verbose) {
              Rcpp::Rcout << "Checking adjacent ID: " << id << " at distance: " << distance << std::endl;
            }
            if (orderMap.at(id) < read_order) {
              if (verbose) {
                Rcpp::Rcout << "Checking id: " << id << std::endl;
              }
              if ((lowerBounds.find(id) != lowerBounds.end() && item.second[0] <= lowerBounds.at(id)) || id == "poly_a" || id == "poly_t" || id == "barcode" || id == "umi") {
                if (!found_end) { // Only update new_start if found_start is false
                  new_end = item.second[1] - 1;
                }
                found_end = true; // Make sure to update the flag
                if(verbose) {
                  Rcpp::Rcout << "New stop found: " << new_end << std::endl;
                }
              }
              continue;
            } else {
              if ((lowerBounds.find(id) != lowerBounds.end() && item.second[0] <= lowerBounds.at(id)) || id == "poly_a" || id == "poly_t" || id == "barcode" || id == "umi") {
                if (!found_start) { // Only update new_start if found_start is false
                  new_start = item.second[2] + 1;
                }
                found_start = true;
                if(verbose) {
                  Rcpp::Rcout << "New start found: " << new_start << std::endl;
                }
              }
              continue;
            }
          }
        }
      }
      if (new_start == -1 && new_end != -1) {
        if(verbose){
          Rcpp::Rcout << "New start not found, setting to sentinel value 0" << std::endl;
        }
        new_start = 1;  // Note the corrected assignment operator
        break;
      }
      if (found_start && found_end) {
        if(verbose) {
          Rcpp::Rcout << "Both new start and end found, exiting loop." << std::endl;
        }
        break;
      }
      if (!found_start || !found_end) {
        search_distance++;
        if(verbose) {
          Rcpp::Rcout << "New position not found, increasing search distance to: " << search_distance << std::endl;
        }
      }
    }
    if (new_start != -1 && new_end != -1) {
      new_read_coords = {0, new_start, new_end};
      if(verbose){
        Rcpp::Rcout << "putative new coordinates are " << new_read_coords[0] <<" "<< new_read_coords[1] << " " << new_read_coords[2] << std::endl;
      }
      parsedSignature[read_key] = new_read_coords;
    }
    std::stringstream newSignature;
    for (const auto& item : parsedSignature) {
      newSignature << "|" << item.first << ":" << item.second[0] << ":" << item.second[1] << ":" << item.second[2];
    }
    newSignature << "<" << category << ">";
    if(verbose){
      Rcpp::Rcout << "putative new signature is " << newSignature.str() <<  "!" << std::endl;
    }
    return newSignature.str().substr(1);}
  // Check for "rescuable_concatenates"
  else if (parsedSignature.at("forw_primer")[0] <= lowerBounds.at("forw_primer") &&
    parsedSignature.at("rc_forw_primer")[0] <= lowerBounds.at("rc_forw_primer") && 
    parsedSignature.find("read") != parsedSignature.end() &&
    parsedSignature.find("rc_read") != parsedSignature.end()) {
    category = "rescuable_concatenates";
    if(verbose){
      Rcpp::Rcout << "Mapped as a rescuable concatenate, now attempting to rescue!" << std::endl;
      Rcpp::Rcout << barcoded_sigstring << std::endl;
    }
    // Create a vector to hold the sorted elements
    std::vector<std::tuple<std::string, int, int, int>> sortedElements;
    for (const auto& item : parsedSignature) {
      std::string id = item.first;
      if (id != "read" && id != "rc_read") {
        sortedElements.push_back(std::make_tuple(id, orderMap.at(id), item.second[1], item.second[2]));
      }
    }
    // Debug: Print elements before sorting
    if (verbose) {
      Rcpp::Rcout << "Elements before sorting:" << std::endl;
      for (const auto& item : sortedElements) {
        Rcpp::Rcout << std::get<0>(item) << " : " << std::get<2>(item) << "-" << std::get<3>(item) << std::endl;
      }
    }
    // Sort by start position
    std::sort(sortedElements.begin(), sortedElements.end(),
      [](const auto& a, const auto& b) {
        return std::get<2>(a) < std::get<2>(b);
      });
    // Debug: Print elements after sorting
    if (verbose) {
      Rcpp::Rcout << "Elements after sorting:" << std::endl;
      for (const auto& item : sortedElements) {
        Rcpp::Rcout << std::get<0>(item) << " : " << std::get<2>(item) << "-" << std::get<3>(item) << std::endl;
      }
    }
    int read_start = parsedSignature.at("read")[2] + 1;
    int read_end = parsedSignature.at("rc_read")[1] - 1;
    if(verbose){
      Rcpp::Rcout << "read start is " << read_start << " and read end is " << read_end << std::endl;
    }
    // Find read start and end based on sorted elements
    for (const auto& item : sortedElements) {
      std::string id = std::get<0>(item);
      int start = std::get<2>(item);
      int end = std::get<3>(item);
      if (start > read_start) {
        read_end = start - 1;
        break;
      }
      if (end < read_end) {
        read_start = end + 1;
      }
    }
    // Debug: Print new read start and end
    if (verbose) {
      Rcpp::Rcout << "New read start: " << read_start << ", New read end: " << read_end << std::endl;
    }
    // Update parsedSignature for read and rc_read
    parsedSignature["read"] = {0, read_start, read_end};
    parsedSignature["rc_read"] = {0, read_end + 1, parsedSignature.at("rc_read")[2]};
    // Step 4: Re-sort the complete list, including updated read and rc_read
    sortedElements.push_back(std::make_tuple("read", orderMap.at("read"), read_start, read_end));
    sortedElements.push_back(std::make_tuple("rc_read", orderMap.at("rc_read"), read_end + 1, parsedSignature.at("rc_read")[2]));
    std::sort(sortedElements.begin(), sortedElements.end(),
      [](const auto& a, const auto& b) {
        return std::get<2>(a) < std::get<2>(b);
      });
    if (verbose) {
      Rcpp::Rcout << "Elements after sorting:" << std::endl;
      for (const auto& item : sortedElements) {
        Rcpp::Rcout << std::get<0>(item) << " : " << std::get<2>(item) << "-" << std::get<3>(item) << std::endl;
      }
    }
    // Initialize variables for the new start and stop positions
    int new_read_start = -1;
    int new_read_end = -1;
    int new_rc_read_start = -1;
    int new_rc_read_end = -1;
    // Iterate through the sorted elements to update start and stop positions
    for (size_t i = 0; i < sortedElements.size(); ++i) {
      if (std::get<0>(sortedElements[i]) == "read") {
        new_read_start = std::get<3>(sortedElements[i - 1]) + 1;  // Stop of the element before 'read' + 1
        new_read_end = std::get<2>(sortedElements[i + 1]) - 1;  // Start of the element after 'read' - 1
      }
      if (std::get<0>(sortedElements[i]) == "rc_read") {
        new_rc_read_start = std::get<3>(sortedElements[i - 1]) + 1;  // Stop of the element before 'rc_read' + 1
        new_rc_read_end = std::get<2>(sortedElements[i + 1]) - 1;  // Start of the element after 'rc_read' - 1
      }
    }
    // Update parsedSignature with new start and stop positions
    parsedSignature["read"] = {0, new_read_start, new_read_end};
    parsedSignature["rc_read"] = {0, new_rc_read_start, new_rc_read_end};
    // Debug: Print new read and rc_read start and end
    if (verbose) {
      Rcpp::Rcout << "New read start: " << new_read_start << ", New read end: " << new_read_end << std::endl;
      Rcpp::Rcout << "New rc_read start: " << new_rc_read_start << ", New rc_read end: " << new_rc_read_end << std::endl;
    }
    // Find and update 'read' and 'rc_read' in sortedElements
    for (auto& item : sortedElements) {
      if (std::get<0>(item) == "read") {
        std::get<2>(item) = new_read_start;
        std::get<3>(item) = new_read_end;
      }
      if (std::get<0>(item) == "rc_read") {
        std::get<2>(item) = new_rc_read_start;
        std::get<3>(item) = new_rc_read_end;
      }
    }
    // Generate the new concatenated signature
    std::stringstream newSignature;
    for (const auto& item : sortedElements) {
      newSignature << "|" << std::get<0>(item) << ":" << 0 << ":" << std::get<2>(item) << ":" << std::get<3>(item);
    }
    // Debug: Print new concatenated signature
    if (verbose) {
      Rcpp::Rcout << "New concatenated signature: " << newSignature.str() << std::endl;
    }
    newSignature << "<" << category << ">";
    return newSignature.str().substr(1);
  }
  // Check for "parallel_f_concatenates"
  if (category == "forwards" && idCount["forw_primer"] > 1) {
    category = "parallel_f_concatenates";
  }
  // Check for "parallel_r_concatenates"
  if (category == "reverses" && idCount["rc_forw_primer"] > 1) {
    category = "parallel_r_concatenates";
  }
  return barcoded_sigstring + "<" + category + ">";
}

std::string read_annotator(const std::string& signatureString, const std::unordered_map<std::string, int>& orderMap) {
  std::vector<std::pair<int, std::string>> orderedComponents;
  // Split the signature string by '|'
  std::vector<std::string> components = split(signatureString, '|');
  for (const auto& component : components) {
    // Extract the ID from each component (the part before the first ':')
    std::string id = component.substr(0, component.find(":"));
    // Check if this ID is in our order map
    if (orderMap.find(id) != orderMap.end()) {
      int order = orderMap.at(id);
      orderedComponents.push_back({order, component});
    }
  }
  // Sort by the order
  std::sort(orderedComponents.begin(), orderedComponents.end());
  std::stringstream reorderedSigString;
  int order_of_read = orderMap.at("read");
  int order_of_rc_read = orderMap.at("rc_read");
  for (const auto& item : orderedComponents) {
    reorderedSigString << "|" << item.second;
    // If the order of the current item is one less than the order of "read", then insert "read"
    if (item.first == order_of_read-1) {
      // Inserting a read with arbitrary start and stop positions (0, 100) for now
      // 
      reorderedSigString << "|read:0:0:100";
    }
    // If the order of the current item is one less than the order of "rc_read", then insert "rc_read"
    if (item.first == order_of_rc_read-1) {
      // Inserting an rc_read with arbitrary start and stop positions (0, 100) for now
      reorderedSigString << "|rc_read:0:0:100";
    }
  }
  return reorderedSigString.str().substr(1);  // Remove the leading "|"
}

std::string bajbatch(const std::string sigstring, 
  const std::unordered_map<std::string, std::unordered_map<std::string, std::tuple<int, int, std::string>>>& layoutMap, 
  const std::unordered_map<std::string, int>& orderMap,
  const std::unordered_map<std::string, double>& lowerBounds, 
  bool verbose = false) {
  std::stringstream barcoded_sigstring;
  std::unordered_map<std::string, std::vector<int>> parsed_sigstring = sig_parse(sigstring);
  for (const auto& annot : parsed_sigstring) {
    std::string id = annot.first;
    std::vector<int> coords = annot.second;
    int edit_distance = coords[0];
    int start_pos = coords[1];
    int stop_pos = coords[2];
    if(verbose) {
      Rcpp::Rcout << "Checking piece identity!"<< std::endl;
      Rcpp::Rcout << id << std::endl;
    }
    // Skip this iteration if the edit distance is greater than the lower bound
    if (lowerBounds.find(id) != lowerBounds.end() && edit_distance > lowerBounds.at(id)) {
      barcoded_sigstring << "|" << id << ":" << edit_distance << ":" << start_pos << ":" << stop_pos;
      if(verbose) {
        Rcpp::Rcout << "Skipping this id: " << id << std::endl;
        Rcpp::Rcout << "Updated Signature so far: " << barcoded_sigstring.str() << std::endl;
      }
      continue;
    }
    // Check if this ID is a flanking adapter
    if (layoutMap.find(id) != layoutMap.end()) {
      std::string direction = layoutMap.at(id).begin()->first;
      std::string type = std::get<2>(layoutMap.at(id).begin()->second);
      if (type == "flanking_adapter") {
        std::string barcode_id = (direction == "forward") ? "barcode" : "rc_barcode";
        std::string umi_id = (direction == "forward") ? "umi" : "rc_umi";
        // Annotate barcode
        if (layoutMap.find(barcode_id) != layoutMap.end() && layoutMap.at(barcode_id).find(direction) != layoutMap.at(barcode_id).end()) {
          int barcode_length = std::get<0>(layoutMap.at(barcode_id).at(direction));
          int new_start, new_stop;
          if (direction == "forward") {
            new_start = stop_pos + 1;
            new_stop = new_start + barcode_length - 1;
          }
          if (direction == "reverse") {
            new_stop = start_pos - 1;
            new_start = new_stop - barcode_length + 1;
          }
          barcoded_sigstring << "|" << barcode_id << ":" << 0 << ":" << new_start << ":" << new_stop;
        }
        // Annotate UMI
        if (layoutMap.find(umi_id) != layoutMap.end() && layoutMap.at(umi_id).find(direction) != layoutMap.at(umi_id).end()) {
          int barcode_length = std::get<0>(layoutMap.at(barcode_id).at(direction));
          int umi_length = std::get<0>(layoutMap.at(umi_id).at(direction));
          int new_start, new_stop;
          
          if (direction == "forward") {
            new_start = stop_pos + barcode_length + 1;
            new_stop = new_start + umi_length - 1;
          } else {
            new_stop = start_pos - barcode_length - 1;
            new_start = new_stop - umi_length + 1;
          }
          barcoded_sigstring << "|" << umi_id << ":" << 0 << ":" << new_start << ":" << new_stop;
        }
      } 
    }
    // Append the existing id and coordinates
    barcoded_sigstring << "|" << id << ":" << edit_distance << ":" << start_pos << ":" << stop_pos;
    if(verbose){
      Rcpp::Rcout << "Updated Signature with barcodes and umis:  " << barcoded_sigstring.str().substr(1) << std::endl;
    }
  }
  std::string ordered_sigstring = read_annotator(barcoded_sigstring.str().substr(1), orderMap);
  if(verbose){
    std::string category = error_categorization(ordered_sigstring, lowerBounds,layoutMap, orderMap, true);
    return category;
  } else {
    std::string category = error_categorization(ordered_sigstring, lowerBounds,layoutMap, orderMap, false);
    return category;
  }
}

// [[Rcpp::export]]
Rcpp::StringVector revcompR(Rcpp::StringVector sequences) {
  int n = sequences.size();
  Rcpp::StringVector results(n);
  int j = 0; // Counter for results
  for (int i = 0; i < n; ++i) {
    if (sequences[i] == NA_STRING) {
      results[j] = NA_STRING;
      continue;
    }
    results[j] = revcomp(Rcpp::as<std::string>(sequences[i]));
    ++j;
  }
  return results;
}

// [[Rcpp::export]]
LogicalVector check_whitelist(DataFrame r_whitelist_df, NumericVector mutations) {
  NumericVector whitelist_bcs = r_whitelist_df["whitelist_bcs"];
  std::unordered_set<int64_t> whitelist_set;
  for (int i = 0; i < whitelist_bcs.size(); ++i) {
    int64_t whitelist_entry = static_cast<int64_t>(whitelist_bcs[i]);
    whitelist_set.insert(whitelist_entry);
  }
  LogicalVector in_whitelist(mutations.size());
  for (int i = 0; i < mutations.size(); ++i) {
    int64_t mutated_barcode = static_cast<int64_t>(mutations[i]);
    in_whitelist[i] = (whitelist_set.find(mutated_barcode) != whitelist_set.end());
  }
  return in_whitelist;
}

// [[Rcpp::export]]
Rcpp::NumericVector barcodes_to_bits(Rcpp::StringVector barcodes) {
  int n = barcodes.size();
  std::vector<int64_t> results(n);  // Changed to int64_t
  for(int i = 0; i < n; ++i) {
    std::string barcode = Rcpp::as<std::string>(barcodes[i]);
    int64_t result = 0;  // Changed to int64_t
    for (char &c : barcode) {
      result <<= 2;
      switch (c) {
      case 'A': break;
      case 'C': result |= 1; break; 
      case 'T': result |= 2; break;
      case 'G': result |= 3; break;
      }
    } 
    results[i] = result;
  } 
  return Rcpp::toInteger64(results);
}

// [[Rcpp::export]]
Rcpp::StringVector bits_to_barcodes(Rcpp::NumericVector input, int barcode_length = 16, bool verbose = false) {
  // Convert the NumericVector to std::vector<int64_t>
  std::vector<int64_t> cpp_bits = Rcpp::fromInteger64(input);
  int n = cpp_bits.size();
  Rcpp::StringVector results(n);
  for(int i = 0; i < n; ++i) {
    int64_t int64_code = cpp_bits[i];  // Changed to int64_t and used cpp_bits
    std::string result;
    for (int j = 0; j < barcode_length; ++j) {
      if (verbose) {
        Rcpp::Rcout << "Current int64_code: " << int64_code << std::endl;
        Rcpp::Rcout << "Last 2 bits: " << (int64_code & 3) << std::endl;
      }
      switch (int64_code & 3) {
      case 0: result.insert(result.begin(), 'A'); break;
      case 1: result.insert(result.begin(), 'C'); break;
      case 2: result.insert(result.begin(), 'T'); break;
      case 3: result.insert(result.begin(), 'G'); break;
      }
      int64_code >>= 2;
    }
    results[i] = result;
  }
  return results;
}

// [[Rcpp::export]]
Rcpp::NumericVector mutate_and_check(Rcpp::NumericVector barcode, int max_mutations, Rcpp::DataFrame r_whitelist_df, int barcode_length = 16) {
  // Convert R integer64 barcode to int64_t using RcppInt64
  std::vector<int64_t> barcodes_vec = Rcpp::fromInteger64(barcode);
  int64_t int64_barcode = barcodes_vec[0];  // Assuming you're working with a single barcode
  Rcpp::NumericVector r_whitelist = r_whitelist_df["whitelist_bcs"];
  std::unordered_set<int64_t> whitelist;
  std::vector<int64_t> whitelist_cpp = Rcpp::fromInteger64(r_whitelist);
  for (const auto& val : whitelist_cpp) {
    whitelist.insert(val);
  }
  int num_bits = 2 * barcode_length;
  std::set<int64_t> unique_results;
  recursive_flip_bits_list(int64_barcode, int64_barcode, 0, max_mutations, num_bits, unique_results, whitelist);
  // Convert the set back to a NumericVector for R using RcppInt64
  Rcpp::NumericVector results = Rcpp::toInteger64(std::vector<int64_t>(unique_results.begin(), unique_results.end()));
  return results;
}

// [[Rcpp::export]]
Rcpp::DataFrame baj_extract(std::vector<std::string>& sigstrings,Rcpp::DataFrame whitelist_df, Rcpp::DataFrame df, 
  bool verbose = false, int max_iterations = 4, int max_mutations = 4, 
  int barcode_length = 16, bool barcorrect = false, int nthreads = 1, bool jaccard_on = true) {
  Rcpp::StringVector ids = df["id"];
  Rcpp::StringVector fastq_files = df["fastq_files"];
  Rcpp::StringVector qcs = df["qc"];
  std::vector<std::string> extracted_ids, extracted_barcodes, extracted_bcqcs, extracted_umis, extracted_reads, extracted_qcs;
  std::vector<std::string> updated_sig_ids;
  // Create a set for the whitelisted barcodes in bitwise representation
  Rcpp::NumericVector r_whitelist = whitelist_df["whitelist_bcs"];
  std::unordered_set<int64_t> whitelist_set;
  std::vector<int64_t> whitelist_cpp = Rcpp::fromInteger64(r_whitelist);
  for (const auto& val : whitelist_cpp) {
    whitelist_set.insert(val);
  }
  // Iterate over each signature string
  for (size_t i = 0; i < sigstrings.size(); ++i) {
    // Parse the signature string
    std::unordered_map<std::string, std::vector<int>> parsedSignature = sig_parse(sigstrings[i]);
    size_t pos = sigstrings[i].find_last_of("<");
    std::string category = sigstrings[i].substr(pos);
    // Function to extract the subsequence based on parsedSignature
    auto extract_subsequence = [&](std::string id_name, std::string sequence) -> std::string {
      int start = parsedSignature[id_name][1] - 1;
      int end = parsedSignature[id_name][2];
      if (start == 0){
        start = 1;
      }
      if(end > sequence.size()){
        return category = "undecided";
      }
      if((end-start) < 5){
        return category = "undecided";
      }
      if (start >= sequence.size() || end > sequence.size()) {
        return category = "undecided";
      }
      return sequence.substr(start, end - start);
    };
    sigstrings[i] = Rcpp::as<std::string>(ids[i]) + "*" + sigstrings[i];
    if(verbose){
      Rcpp::Rcout << sigstrings[i] << std::endl;
    }
    if(category == "undecided"||category == "<undecided>"){
      continue;
    }
    if (category == "<forward>" || category == "<reverse>") {
      std::string prefix = (category == "<forward>") ? "" : "rc_";
      int start_pos = parsedSignature[prefix + "read"][1];
      int stop_pos = parsedSignature[prefix + "read"][2];
      // Set stop_pos to the length of the fastq_file string if stop_pos is zero
      if (category == "<forward>" && stop_pos == 0) {
        stop_pos = fastq_files[i].size();
      }
      // Skip this iteration if start_pos > stop_pos
      if (start_pos > stop_pos) {
        if(verbose) {
          Rcpp::Rcout << "Skipping: start_pos > stop_pos for ID " << ids[i] << std::endl;
        }
        continue;
      }
      std::string read = Rcpp::as<std::string>(fastq_files[i]).substr(start_pos, stop_pos - start_pos);
      std::string qc = Rcpp::as<std::string>(qcs[i]).substr(start_pos, stop_pos - start_pos);
      extracted_ids.push_back((sigstrings[i]) + category);
      if(category == "<reverse>"){
        std::string rc_barcode = extract_subsequence("rc_barcode", Rcpp::as<std::string>(fastq_files[i]));
        std::string revcomp_barcode = revcomp(rc_barcode);
        extracted_barcodes.push_back(revcomp_barcode);
        std::string rc_umi = extract_subsequence("rc_umi", Rcpp::as<std::string>(fastq_files[i]));
        std::string revcomp_umi = revcomp(rc_umi);
        extracted_umis.push_back(revcomp_umi);
        std::string bcqc = extract_subsequence("rc_barcode", Rcpp::as<std::string>(qcs[i]));
        extracted_bcqcs.push_back(bcqc);
      } else {
        extracted_barcodes.push_back(extract_subsequence(prefix + "barcode", Rcpp::as<std::string>(fastq_files[i])));
        extracted_umis.push_back(extract_subsequence(prefix + "umi", Rcpp::as<std::string>(fastq_files[i])));
        extracted_bcqcs.push_back(extract_subsequence(prefix + "barcode", Rcpp::as<std::string>(qcs[i])));
      }
      extracted_reads.push_back(read);
      extracted_qcs.push_back(qc);
    }
    if (category == "<rescuable_concatenates>") {
      // First split
      std::string extracted_read = extract_subsequence("read", Rcpp::as<std::string>(fastq_files[i]));
      if(extracted_read.size() < 30){
        continue;
      }
      extracted_ids.push_back(sigstrings[i] + "_splitid=<forward>");
      extracted_barcodes.push_back(extract_subsequence("barcode", Rcpp::as<std::string>(fastq_files[i])));
      extracted_umis.push_back(extract_subsequence("umi", Rcpp::as<std::string>(fastq_files[i])));
      extracted_reads.push_back(extract_subsequence("read", Rcpp::as<std::string>(fastq_files[i])));
      extracted_qcs.push_back(extract_subsequence("read", Rcpp::as<std::string>(qcs[i])));
      extracted_bcqcs.push_back(extract_subsequence("barcode", Rcpp::as<std::string>(qcs[i])));
      
      // Second split
      std::string extracted_rc_read = extract_subsequence("rc_read", Rcpp::as<std::string>(fastq_files[i]));
      if(extracted_rc_read.size() < 30){
        continue;
      }
      extracted_ids.push_back(sigstrings[i] + "_splitid=<reverse>");
      std::string rc_barcode = extract_subsequence("rc_barcode", Rcpp::as<std::string>(fastq_files[i]));
      std::string revcomp_barcode = revcomp(rc_barcode);
      extracted_barcodes.push_back(revcomp_barcode);
      std::string rc_umi = extract_subsequence("rc_umi", Rcpp::as<std::string>(fastq_files[i]));
      std::string revcomp_umi = revcomp(rc_umi);
      extracted_umis.push_back(revcomp_umi);
      extracted_reads.push_back(extract_subsequence("rc_read", Rcpp::as<std::string>(fastq_files[i])));
      extracted_qcs.push_back(extract_subsequence("rc_read", Rcpp::as<std::string>(qcs[i])));
      std::string bcqc = extract_subsequence("rc_barcode", Rcpp::as<std::string>(qcs[i]));
      extracted_bcqcs.push_back(bcqc);
    }
  }
  
  std::vector<std::string> wl_barcodes, nwl_barcodes;
  std::vector<std::string> wl_umis, nwl_umis;
  std::vector<std::string> wl_reads, nwl_reads;
  std::vector<std::string> wl_qcs, nwl_qcs;
  std::vector<std::string> wl_ids, nwl_ids;
  std::vector<std::string> wl_bcqcs, nwl_bcqcs;
  
  std::vector<int64_t> bit_barcodes = barcodes_to_bits_cpp(extracted_barcodes);
  for (size_t i = 0; i < extracted_barcodes.size(); ++i) {
    int64_t bit_barcode = bit_barcodes[i];
    if (whitelist_set.find(bit_barcode) != whitelist_set.end()) {
      wl_barcodes.push_back(extracted_barcodes[i]);
      wl_umis.push_back(extracted_umis[i]);
      wl_reads.push_back(extracted_reads[i]);
      wl_qcs.push_back(extracted_qcs[i]);
      wl_ids.push_back(extracted_ids[i]);
      wl_bcqcs.push_back(extracted_bcqcs[i]);
    } else {
      nwl_barcodes.push_back(extracted_barcodes[i]);
      nwl_umis.push_back(extracted_umis[i]);
      nwl_reads.push_back(extracted_reads[i]);
      nwl_qcs.push_back(extracted_qcs[i]);
      nwl_ids.push_back(extracted_ids[i]);
      nwl_bcqcs.push_back(extracted_bcqcs[i]);
    }
  }
  // Perform barcode error correction only if barcorrect is true
  if (barcorrect){
    omp_set_num_threads(nthreads);
    std::vector<int64_t> nwl_intcodes = barcodes_to_bits_cpp(nwl_barcodes);
    // OpenMP parallelization starts here
#pragma omp parallel for
    for (size_t i = 0; i < nwl_barcodes.size(); ++i) {
      int64_t original_bit_barcode = nwl_intcodes[i];
      int min_hamming_distance = max_mutations;
      std::set<int64_t> min_hamming_results;
      // Ensure this function is thread-safe
      recursive_flip_bits_cpp(original_bit_barcode, original_bit_barcode, 0, max_iterations, 2 * barcode_length, whitelist_set, 
        min_hamming_distance, min_hamming_results);
      // If a single barcode is found with the minimum Hamming distance, correct the original
      if(verbose){
#pragma omp critical
{
  Rcpp::Rcout << "Original barcode: " << nwl_barcodes[i] << " Number of hamming-distance mapped barcodes found in whitelist: " << min_hamming_results.size() << std::endl;
  if(min_hamming_results.size() <= 5){
    std::vector<int64_t> min_hamming_results_vec(min_hamming_results.begin(), min_hamming_results.end());
    std::vector<std::string> results = bits_to_barcodes_cpp(min_hamming_results_vec);
    Rcpp::Rcout << "Generated barcodes: ";
    for (const auto& str : results){
      Rcpp::Rcout << str << ' ';
    }
    Rcpp::Rcout << std::endl;
    // Iterate over each generated barcode to find mismatch positions
    for (const auto& generated_barcode : results){
      std::vector<int> mismatch_positions = find_mismatch_positions(nwl_barcodes[i], generated_barcode);
      Rcpp::Rcout << "Mismatch positions for " << generated_barcode << ": ";
      for (const auto& pos : mismatch_positions){
        Rcpp::Rcout << pos << ' ';
      }
      Rcpp::Rcout << std::endl;
    }
  }
}
      }
      if (min_hamming_results.size() == 1) {
        int64_t corrected_bit_barcode = *min_hamming_results.begin();
        std::string corrected_barcode = bits_to_barcodes_cpp({corrected_bit_barcode})[0];
        nwl_barcodes[i] = corrected_barcode;
        nwl_ids[i] += "{orig_hamming_" + nwl_barcodes[i] + "_" + corrected_barcode + ":" + std::to_string(min_hamming_distance) + "}";
        if(verbose){
#pragma omp critical
{
  Rcpp::Rcout << "Hamming corrected barcode " << nwl_barcodes[i] << " to " << corrected_barcode << std::endl;
}
        }
      }
      if (min_hamming_results.size() <= 20 & jaccard_on) {
        std::vector<int64_t> min_hamming_results_vec(min_hamming_results.begin(), min_hamming_results.end());
        std::vector<std::string> candidate_barcodes = bits_to_barcodes_cpp(min_hamming_results_vec);
        auto result = best_matches_jaccard_distance_qgram(nwl_barcodes[i], candidate_barcodes, 2);
        std::vector<std::string> best_matches = result.first;
        double min_jaccard_distance = result.second;
        if (best_matches.size() == 1) {
          nwl_barcodes[i] = best_matches[0];
          nwl_ids[i] += "{orig_jaccard_" + nwl_barcodes[i] + "_" + best_matches[0] + ":" + std::to_string(min_jaccard_distance) + "}";
        }
        if (verbose) {
#pragma omp critical
{
  Rcpp::Rcout << "Jaccard corrected barcode " << nwl_barcodes[i] << " to ";
  for (const auto& match : best_matches) {
    Rcpp::Rcout << match << ' ';
  }
  Rcpp::Rcout << std::endl;
}
        }
      }
    }
    wl_barcodes.insert(wl_barcodes.end(), nwl_barcodes.begin(), nwl_barcodes.end());
    wl_umis.insert(wl_umis.end(), nwl_umis.begin(), nwl_umis.end());
    wl_reads.insert(wl_reads.end(), nwl_reads.begin(), nwl_reads.end());
    wl_qcs.insert(wl_qcs.end(), nwl_qcs.begin(), nwl_qcs.end());
    wl_ids.insert(wl_ids.end(), nwl_ids.begin(), nwl_ids.end());
    wl_bcqcs.insert(wl_bcqcs.end(), nwl_bcqcs.begin(), nwl_bcqcs.end());
  }
  return Rcpp::DataFrame::create(Rcpp::Named("sig_id") = wl_ids,
    Rcpp::Named("barcode") = wl_barcodes,
    Rcpp::Named("barcode_qc") = wl_bcqcs,
    Rcpp::Named("umi") = wl_umis,
    Rcpp::Named("filtered_read") = wl_reads,
    Rcpp::Named("filtered_qc") = wl_qcs);
}

// [[Rcpp::export]]
StringVector bajbatch(DataFrame null_distance, DataFrame read_layout, StringVector sigstrings, bool verbose = false) {
  if (verbose) {
    Rcpp::Rcout << "Null Distance DataFrame: " << std::endl;
    CharacterVector query_id = null_distance["query_id"];
    NumericVector null_dist = null_distance["null_distance"];
    NumericVector sd_null = null_distance["sd_null"];
    for (int i = 0; i < query_id.size(); ++i) {
      Rcpp::Rcout << as<std::string>(query_id[i]) << ", " << null_dist[i] << ", " << sd_null[i] << std::endl;
    }
  }
  // Calculate the lower bounds for each id
  std::unordered_map<std::string, double> lowerBounds = calculateLowerBounds(null_distance);
  if(verbose) {
    Rcpp::Rcout << "lowerBounds: " << std::endl;
    for (const auto& pair : lowerBounds) {
      Rcpp::Rcout << "ID: " << pair.first 
                  << ", Value: " << pair.second << std::endl;
    }
  }
  // Convert read_layout to an appropriate C++ data structure
  std::unordered_map<std::string, std::unordered_map<std::string, 
    std::tuple<int, int, std::string>>> layoutMap = convertReadLayout(read_layout);  
  std::unordered_map<std::string, int> orderMap = generateOrderMap(read_layout);
  if (verbose) {
    Rcpp::Rcout << "Layout Map: " << std::endl;
    for (const auto& outer_pair : layoutMap) {
      for (const auto& inner_pair : outer_pair.second) {
        Rcpp::Rcout << "ID: " << outer_pair.first 
                    << ", Direction: " << inner_pair.first
                    << ", Expected Length: " << std::get<0>(inner_pair.second)
                    << ", Order: " << std::get<1>(inner_pair.second) << std::endl
                    << ", Type: " << std::get<2>(inner_pair.second) << std::endl;
      }
    }
  }
  // Iterate through each signature string to parse and update it
  for (int i = 0; i < sigstrings.size(); ++i) {
    std::string sigstring = as<std::string>(sigstrings[i]);
    if(verbose){
      std::string output = bajbatch(sigstring, layoutMap, orderMap, lowerBounds, true);
      sigstrings[i] = output;
    } else {
      std::string output = bajbatch(sigstring, layoutMap, orderMap, lowerBounds, false);
      sigstrings[i] = output;
    }
  }
  return sigstrings;
}
