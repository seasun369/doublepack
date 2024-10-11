#include <catch2/catch.hpp>
#include <iostream>

#include "../PSS/pss.h"
#include "../RMFE/RMFE.h"

using namespace std;
using namespace packed_shamir;

void run_pss_test(int members, int packed_number, int poly_degree, int corrupt_parties) {
    std::cout << "\n==== Test Parameters ====" << std::endl;
    std::cout << "Members: " << members << std::endl;
    std::cout << "Packed number: " << packed_number << std::endl;
    std::cout << "Polynomial degree: " << poly_degree << std::endl;
    std::cout << "Corrupt parties: " << corrupt_parties << std::endl;
    std::cout << "=========================" << std::endl;

    ZZ gr_prime = conv<ZZ>(2);  //p 
    long positive_integer = 4; //k
    long gr_degree = 4;         //r

    gr galoisring(gr_prime, positive_integer, gr_degree);

    // Check the polynomial degree condition before creating the scheme
    if (poly_degree < packed_number + corrupt_parties - 1 || poly_degree > members - packed_number) {
        throw std::invalid_argument("Invalid polynomial degree");
    }

    packed_shamir::scheme Scheme(members, packed_number, poly_degree, galoisring);

    vec_ZZ_pE test;
    test.SetLength(packed_number);

    for (int i = 0; i < packed_number; i++) {
        test[i] = random_ZZ_pE();
    }

    vec_ZZ_pE v = Scheme.create_shares(test);
    vector<int> party;
    for (int i = 1; i <= members; i++) {
        party.push_back(i);
    }

    vec_ZZ_pE u = Scheme.packed_reconstruct_shares(party, v);

    for (int i = 0; i < packed_number; i++) {
        REQUIRE(test[i] == u[i]);
    }

    std::cout << "Test completed successfully." << std::endl;
}

TEST_CASE("PSS with minimum valid polynomial degree", "[pss][degrees][min]") {
    int members = 10;
    int packed_number = 3;
    int corrupt_parties = 4;
    int min_poly_degree = packed_number + corrupt_parties - 1;
    REQUIRE_NOTHROW(run_pss_test(members, packed_number, min_poly_degree, corrupt_parties));
}

TEST_CASE("PSS with maximum valid polynomial degree", "[pss][degrees][max]") {
    int members = 10;
    int packed_number = 3;
    int corrupt_parties = 4;
    int max_poly_degree = members - packed_number;
    REQUIRE_NOTHROW(run_pss_test(members, packed_number, max_poly_degree, corrupt_parties));
}

TEST_CASE("PSS with middle valid polynomial degree", "[pss][degrees][mid]") {
    int members = 10;
    int packed_number = 3;
    int corrupt_parties = 4;
    int mid_poly_degree = (packed_number + corrupt_parties - 1 + members - packed_number) / 2;
    REQUIRE_NOTHROW(run_pss_test(members, packed_number, mid_poly_degree, corrupt_parties));
}

TEST_CASE("PSS with invalid polynomial degree - too low", "[pss][degrees][invalid][low]") {
    int members = 10;
    int packed_number = 3;
    int corrupt_parties = 4;
    int invalid_low_degree = packed_number + corrupt_parties - 2;
    REQUIRE_THROWS_AS(run_pss_test(members, packed_number, invalid_low_degree, corrupt_parties), std::invalid_argument);
}

TEST_CASE("PSS with invalid polynomial degree - too high", "[pss][degrees][invalid][high]") {
    int members = 10;
    int packed_number = 3;
    int corrupt_parties = 4;
    int invalid_high_degree = members - packed_number + 1;
    REQUIRE_THROWS_AS(run_pss_test(members, packed_number, invalid_high_degree, corrupt_parties), std::invalid_argument);
}

// Keep the original test case
TEST_CASE("Original PSS test", "[pss][original]") {
    SECTION("pss"){
        std::cout << "\n==== Original PSS Test ====" << std::endl;
        ZZ gr_prime = conv<ZZ>(2);  //p 
        long positive_integer = 4; //k
        long gr_degree = 4;         //r

        gr galoisring(gr_prime, positive_integer, gr_degree);

        int members = 5;
        int packed_number = 2;
        int poly_degree = 3;
        int corrupt_parties = 2;

        std::cout << "Members: " << members << std::endl;
        std::cout << "Packed number: " << packed_number << std::endl;
        std::cout << "Polynomial degree: " << poly_degree << std::endl;
        std::cout << "Corrupt parties: " << corrupt_parties << std::endl;

        packed_shamir::scheme Scheme(members,packed_number,poly_degree,galoisring);

        vec_ZZ_pE test;
        test.SetLength(2);

        test[0] = random_ZZ_pE();
        test[1] = random_ZZ_pE();

        vec_ZZ_pE v = Scheme.create_shares(test);

        vector<int> party = {1,2,3,4,5};

        vec_ZZ_pE u = Scheme.packed_reconstruct_shares(party, v);

        for(int i=0; i< packed_number; i++)
        {
            REQUIRE(test[i] == u[i]);
        }
        std::cout << "Original PSS test completed successfully." << std::endl;
    }
}