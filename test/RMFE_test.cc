#include <catch2/catch.hpp>
#include <iostream>

#include"../RMFE/RMFE.hpp"
#include"../RMFE/utils.hpp"

using namespace std;

TEST_CASE("RMFE"){
    ZZ p (ZZ(2));
    long k = 16;
    long s = 1;
    long D = 2;
    long n1 = 2;
    long n2 = 2;
    long r1 = 5;
    long r2 = 20;

    vector<long> Input = {9,2,4,1};
    RMFE_GR a(p,s,k,D,n1,n2,r1,r2);

    cout << "Current modulus: " << ZZ_p::modulus() << endl;
    cout << "Current modulus polynomial: " << ZZ_pE::modulus() << endl;

    SECTION("RMFE"){
    a.set_input(Input);
    print(Input);

    a.RMFE_GR_PHI();
    vector<long> result = a.get_result();

    ZZ_pX g_;
    SetCoeff(g_,20,1);
    SetCoeff(g_,8,1);
    SetCoeff(g_,0,1);

    ZZ_pE::init(g_);

    ZZ_pE b = long2ZZpE(result);

    print(result);

    vector<long> res = a.RMFE_GR_PSI(result);

    print(res);

    vector<long> Input2 = {5,9,8,7};
    a.set_input(Input2);
    a.RMFE_GR_PHI();
    vector<long> result2 = a.get_result();

    //ZZ_pE a3 = result * result2;

    //ZZ_pE b = long2ZZpE(result);
    //ZZ_pE c = long2ZZpE(result2);
    //ZZ_pE d = b*c;

    //vector<long> e;

    //ZzpE2Veclong(d,e,r1);

    //vector<long> res3 = a.RMFE_GR_PSI(e);
    //print(res3);
    }

}

