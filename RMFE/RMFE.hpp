#ifndef RMFE_HPP
#define RMFE_HPP


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
#include <iostream>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/ZZ_pE.h>
#include <vector>
#include <fstream>

using namespace std;
using namespace NTL;


#include "utils.hpp"
#include "HenselLift.hpp"
#include "Inverse.hpp"
#include "PrimitiveElement.hpp"


class RMFE_GR{
public:

    ZZ p;
    long s;
    long k;
    long D;
    long n1;
    long n2;
    long r1;
    long r2;
    vector<long> PrimitiveIrred1;
    vector<long> PrimitiveIrred2;


    vector<long> input;
    vector<long> output;

    vector<long> Interpolation1;
    vector<long> Interpolation2;

    
    RMFE_GR() = default;
    RMFE_GR(ZZ p_, long s_, long k_, long D_, long n1_, long n2_,
            long r1_, long r2_);
    
    void set_input(const vector<long>& Input);

    void RMFE_GR_INIT1();
    void RMFE_GR_INIT2();
    void RMFE_GR_INIT2_cache();

    void RMFE_GR_PHI();
    vector<long> RMFE_GR_PHI1(vector<long>& input);
    void RMFE_GR_PHI2(vector<long>& input);

    vector<long> RMFE_GR_PSI(vector<long> Input);

    vector<long> get_result(){ return this->output; }

    void get_phi_kernel(vec_ZZ_pE& vec);

};


#endif
