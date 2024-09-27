#ifndef UTILS_HPP
#define UTILS_HPP

#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/vec_ZZ.h>
#include <NTL/vec_ZZ_p.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/vec_vec_ZZ_pE.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/ZZ_pE.h>
#include <vector>
#include <sstream>
using namespace std;
using namespace NTL;

#include "Inverse.hpp"

vector<long> FindFactor(long n);
ZZ_pX long2ZZpX(const vector<long>& coefficients);
ZZ_pE long2ZZpE(const vector<long>& coefficients);
bool isPowerOfTwo(long n);
void VeczzpE2Veclong(const vec_ZZ_pE& v1, vector<long>& m, long s);
string Veclong2String(const std::vector<long>& vec);
void ZzpE2Veclong(const ZZ_pE& F, vector<long>& m, long s);
void VeczzpE2Vecstring(const vec_ZZ_pE& v1, vector<string>& m, long s);
void ZZpX2long(const ZZ_pX& F, vector<long>& Irred);
void ZZpEX2long(const ZZ_pEX& v1, vector<long>& m, long s);
bool allNonZero(const vec_ZZ_pE& vec);
long nextPowerOf2(long n);
void print(vector<long>& vec);
vector<long> PadVectorToLength(const vector<long>& input, size_t targetLength);
vector<long> TrimVector(const std::vector<long>& input);
vector<long> SplitAndPadVector(const vector<long>& input, long segmentLength, long numSegments);
void Long2ZZpEX(const vector<long>& result, ZZ_pEX& V, long s);
void Long2ZZpEX2(const vector<long>& result, ZZ_pEX& V, long s,long n);
int nearestPerfectSquare(int num);
void fillIrred(const ZZ_pX& F, vector<long>& Irred);
void fillInterpolation(const vec_ZZ_pE& v1, vector<long>& Interpolation, long s);
void FindPrimitivePoly(ZZ_pX& g, ZZ p, long n);
void interpolate_for_GR(ZZ_pEX& f, const vec_ZZ_pE& a, const vec_ZZ_pE& b, ZZ p, long l, long s);
vector<vector<long>> splitVector(const vector<long>& input, int groupSize);
void generateMatrix(const ZZ_pE& F, long d, vector<vector<long>>& matrix);
vector<long> multiplyMatrixByVector(const vector<vector<long>>& matrix, const vector<long>& vec);
long fresh(long x, int k);
vector<vector<long>> compression(const ZZ_pX& F);
ZZ_pE ZZpmulZZpE(const ZZ_p& a, const ZZ_pE& b);
vector<vector<ZZ_p>> ZZpEmatrix2ZZpmatrix(vector<vector<ZZ_pE>>& V, long degree, long k);
vector<vector<long>> matrixMultiply(const vector<vector<long>>& A, const vector<vector<long>>& B);
#endif
