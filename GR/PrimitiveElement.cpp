#include "PrimitiveElement.h"

/*  
    find the primitive element over the galois ring GR(p^k,s) ; 
    user must modify the coefficient of the primitive polynomial, which can be found at https://www.oreilly.com/library/view/built-in-test/9780471624639/18_app.html
*/
ZZ_pE FindPrimitiveElement(ZZ p, long k, long s){

    ZZ q1 = power(p,s);
    ZZ q3 = power(p,k);
    ZZ q6 = power(p,k-1);

    long q4;
    conv(q4,q1-ZZ(1));
    long q5;
    conv(q5,q6);

    ZZ_p::init(p);
    ZZ_pX F;

    /////  user need to modify here
    SetCoeff(F,30,1);
    SetCoeff(F,16,1);
    SetCoeff(F,15,1);
    SetCoeff(F,1,1);
    SetCoeff(F,0,1);
    /////

    ZZ_p::init(q3);
    ZZ_pE::init(F);

    ZZ_pX H;
    SetCoeff(H,1,1);
    ZZ_pE b;
    conv(b,H);

    return power(b,q5);

}
