#include <catch2/catch.hpp>
#include <iostream>

#include"../RMFE/RMFE.hpp"
#include"../RMFE/utils.hpp"

using namespace std;

TEST_CASE("RMFE"){
    SECTION("RMFE"){
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

    a.set_input(Input);
    print(Input);

    a.RMFE_GR_PHI();
    vector<long> result = a.get_result();

    print(result);

    vector<long> res = a.RMFE_GR_PSI(result);

    print(res);

    //vector<long> Input2 = {5,9,8,7};
    //a.set_input(Input2);
    //a.RMFE_GR_PHI();
    //vector<long> result2 = a.get_result();

    //vector<long> c = result * result2;

    //vector<long> res3 = a.RMFE_GR_PSI(result2);
    //print(res3);
    }

}

