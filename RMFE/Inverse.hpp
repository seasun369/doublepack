#ifndef INVERSE_HPP
#define INVERSE_HPP

#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ_pE.h>

using namespace std;
using namespace NTL;

ZZ_pE Inv(ZZ_pE a, long s);

ZZ_pE Inv2(ZZ_pE a, ZZ_pX F, ZZ p, long s, long k);

#endif