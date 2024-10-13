#ifndef DP_H
#define DP_H

#include <vector>
#include "RMFE.h"
#include "scl.h"
#include "gr.h"
#include "pss.h"
namespace dp {
  //Galois ring
  ZZ gr_prime = conv<ZZ>(2);  //p 
  long positive_integer = 4; //k
  long gr_degree = 4;         //r should equal to r2

  gr gring(gr_prime, positive_integer, gr_degree);
  ZZ_pX f = gring.Fps_poly;

  //RMFE, n1*n2 = l
  ZZ p (ZZ(2));
  long k = 16;
  long s = 1;
  long D = 2;
  long n1 = 2;
  long n2 = 2; 
  long r1 = 5;
  long r2 = 20; 

  RMFE_GR rmfe(p,s,k,D,n1,n2,r1,r2);

  // need to redefine
  using FF = ZZ_p;
  using Shr = ZZ_pE;
  using Poly = scl::details::EvPolynomial<FF>;
  using Vec = vector<Shr>;

  template<typename T>
  using vec = std::vector<T>;

}

#endif  // DP_H
