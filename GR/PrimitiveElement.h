#ifndef PRIMITIVEELEMENT_HPP
#define PRIMITIVEELEMENT_HPP

#include <NTL/ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>


using namespace std;
using namespace NTL;

ZZ_pE FindPrimitiveElement(ZZ p, long k, long s);

#endif