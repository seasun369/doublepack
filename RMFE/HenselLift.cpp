#include "HenselLift.hpp"


void HenselLift(ZZ_pX& g_, ZZ_pX& h_, ZZ_pX& s_, ZZ_pX& t_, const ZZ_pX& f, 
                const ZZ_pX& g, const ZZ_pX& h, const ZZ_pX& s, const ZZ_pX& t){
    ZZ_pX e = f - g*h;

    ZZ_pX q,r;
    DivRem(q,r,s*e,h);

    g_ = g + t*e + q*g;
    h_ = h + r;

    ZZ_pX b = s*g_ + t*h_ -1;

    ZZ_pX c,d;
    DivRem(c,d,s*b,h_);

    s_ = s - d;
    t_ = t - t*b - c*g_;

    return;
}

/*
    Hensel Lift of g: g|f in Zp[x], g_|f in Zp^(n+1)[x]
*/
void HenselLift(ZZ_pX& g_, const ZZ_pX& f, const ZZ_pX& g, const ZZ p, long n){
    ZZ_p::init(p);

    ZZ_pX h = f / g;
    ZZ_pX d,s,t;
    XGCD(d,s,t,g,h);

    ZZ_pX g_tmp = g;
    ZZ_pX h_tmp = h;
    ZZ_pX f_tmp = f;
    for(int i = 0; i < n; i++){
        ZZ_pX h_, s_, t_;
        ZZ_p::init(power(p,i+2));
        SetCoeff(f_tmp,0,-1);
        HenselLift(g_,h_,s_,t_,f_tmp,g_tmp,h_tmp,s,t);
        g_tmp = g_;
        h_tmp = h_;
        s = s_;
        t = t_;
    }
    return;

}