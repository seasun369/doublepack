#ifndef GaloisRing
#define GaloisRing

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
#include "PrimitiveElement.h"
using namespace std;
using namespace NTL;

// this file gives the definition of Galois Ring

class gr{
public:
    ZZ p;
    long k;
    long r;
    ZZ_pX Fps_poly;
    vec_ZZ_pE set_T;
    //ZZ_pE element;
    gr(){};
    gr(ZZ p_, long k_, long r_);

    void gr_init_basic_var();
    void get_set_T(vec_ZZ_pE& set_T);
};



#endif