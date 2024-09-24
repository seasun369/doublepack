#ifndef HENSELLIFT_HPP
#define HENSELLIFT_HPP

#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
using namespace std;
using namespace NTL;



void HenselLift(ZZ_pX& g_, const ZZ_pX& f, const ZZ_pX& g, const ZZ p, long n);



#endif