#ifndef DP_H
#define DP_H

#include <vector>
#include "RMFE.h"
#include "scl.h"
#include "gr.h"
#include "pss.h"
namespace dp {

  //RMFE, n1*n2 = l
  inline ZZ p (ZZ(2));
  inline long k = 16;
  inline long s = 1;
  inline long D = 2;
  inline long n1 = 2;
  inline long n2 = 2; 
  inline long r1 = 5;
  inline long r2 = 20; 

  inline RMFE_GR rmfe(p,s,k,D,n1,n2,r1,r2);

    //Galois ring
  inline ZZ gr_prime = conv<ZZ>(2);  //p 
  inline long positive_integer = 16; //k
  inline long gr_degree = 20;         //r should equal to r2


  inline gr gring(gr_prime, positive_integer, gr_degree);
  inline ZZ_pX f = gring.Fps_poly;

  using FF = ZZ_p;
  using Shr = ZZ_pE;
  using Poly = scl::details::EvPolynomial<FF>;
  using Vec = vector<Shr>;

  template<typename T>
  using vec = std::vector<T>;

}

#endif  // DP_H
