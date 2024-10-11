#include <catch2/catch.hpp>
#include <iostream>

#include "../PSS/pss.h"
#include "../RMFE/RMFE.h"

using namespace std;
using namespace packed_shamir;

void print_separator() {
    std::cout << "\n========================================\n" << std::endl;
}

void run_pss_test(int members, int packed_number, int poly_degree, int corrupt_parties) {
    std::cout << "Test Parameters:" << std::endl;
    std::cout << "Members: " << members << ", Packed number: " << packed_number << std::endl;
    std::cout << "Polynomial degree: " << poly_degree << ", Corrupt parties: " << corrupt_parties << std::endl;

    ZZ gr_prime = conv<ZZ>(2);
    long positive_integer = 4;
    long gr_degree = 4;

    gr galoisring(gr_prime, positive_integer, gr_degree);

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
    print_separator();
    std::cout << "==== PSS with minimum valid polynomial degree ====" << std::endl;
    int members = 10, packed_number = 3, corrupt_parties = 4;
    int min_poly_degree = packed_number + corrupt_parties - 1;
    REQUIRE_NOTHROW(run_pss_test(members, packed_number, min_poly_degree, corrupt_parties));
}

TEST_CASE("PSS with maximum valid polynomial degree", "[pss][degrees][max]") {
    print_separator();
    std::cout << "==== PSS with maximum valid polynomial degree ====" << std::endl;
    int members = 10, packed_number = 3, corrupt_parties = 4;
    int max_poly_degree = members - packed_number;
    REQUIRE_NOTHROW(run_pss_test(members, packed_number, max_poly_degree, corrupt_parties));
}

TEST_CASE("PSS with middle valid polynomial degree", "[pss][degrees][mid]") {
    print_separator();
    std::cout << "==== PSS with middle valid polynomial degree ====" << std::endl;
    int members = 10, packed_number = 3, corrupt_parties = 4;
    int mid_poly_degree = (packed_number + corrupt_parties - 1 + members - packed_number) / 2;
    REQUIRE_NOTHROW(run_pss_test(members, packed_number, mid_poly_degree, corrupt_parties));
}

TEST_CASE("PSS with invalid polynomial degree - too low", "[pss][degrees][invalid][low]") {
    print_separator();
    std::cout << "==== PSS with invalid polynomial degree - too low ====" << std::endl;
    int members = 10, packed_number = 3, corrupt_parties = 4;
    int invalid_low_degree = packed_number + corrupt_parties - 2;
    REQUIRE_THROWS_AS(run_pss_test(members, packed_number, invalid_low_degree, corrupt_parties), std::invalid_argument);
}

TEST_CASE("PSS with invalid polynomial degree - too high", "[pss][degrees][invalid][high]") {
    print_separator();
    std::cout << "==== PSS with invalid polynomial degree - too high ====" << std::endl;
    int members = 10, packed_number = 3, corrupt_parties = 4;
    int invalid_high_degree = members - packed_number + 1;
    REQUIRE_THROWS_AS(run_pss_test(members, packed_number, invalid_high_degree, corrupt_parties), std::invalid_argument);
}

TEST_CASE("Original PSS test", "[pss][original]") {
    print_separator();
    std::cout << "==== Original PSS Test ====" << std::endl;
    ZZ gr_prime = conv<ZZ>(2);
    long positive_integer = 4;
    long gr_degree = 4;

    gr galoisring(gr_prime, positive_integer, gr_degree);

    int members = 5, packed_number = 2, poly_degree = 3, corrupt_parties = 2;

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

TEST_CASE("Test packed_create_shares and packed_reconstruct_shares", "[pss][packed]") {
    print_separator();
    std::cout << "==== Test packed_create_shares and packed_reconstruct_shares ====" << std::endl;
    ZZ gr_prime = conv<ZZ>(2);
    long positive_integer = 4;
    long gr_degree = 4;

    gr galoisring(gr_prime, positive_integer, gr_degree);

    int members = 5, packed_number = 2, poly_degree = 3;

    packed_shamir::scheme Scheme(members, packed_number, poly_degree, galoisring);

    vec_ZZ_pE secret;
    secret.SetLength(4);
    for (int i = 0; i < 4; i++) {
        secret[i] = random_ZZ_pE();
    }

    vector<vec_ZZ_pE> shares = Scheme.packed_create_shares(secret);

    REQUIRE(shares.size() == packed_number);

    for (int i = 0; i < packed_number; i++) {
        vector<int> party = {1, 2, 3, 4, 5};
        vec_ZZ_pE reconstructed = Scheme.packed_reconstruct_shares(party, shares[i]);

        REQUIRE(reconstructed.length() == packed_number);
        for (int j = 0; j < packed_number; j++) {
            REQUIRE(reconstructed[j] == secret[i * packed_number + j]);
        }
    }
    std::cout << "Test completed successfully." << std::endl;
}

TEST_CASE("Test create_one_shares and reconstruct_one_shares", "[pss][one_share]") {
    print_separator();
    std::cout << "==== Test create_one_shares and reconstruct_one_shares ====" << std::endl;
    
    ZZ gr_prime = conv<ZZ>(2);
    long positive_integer = 4;
    long gr_degree = 4;

    gr galoisring(gr_prime, positive_integer, gr_degree);

    int members = 5, packed_number = 2, poly_degree = 3;

    packed_shamir::scheme Scheme(members, packed_number, poly_degree, galoisring);

    ZZ_pE secret = random_ZZ_pE();
    long index = 1;

    vec_ZZ_pE shares = Scheme.create_one_shares(secret, index);
    REQUIRE(shares.length() == members);

    ZZ_pE reconstructed = Scheme.reconstruct_one_shares(shares, index);
    REQUIRE(reconstructed == secret);

    std::cout << "Test completed successfully." << std::endl;
}

TEST_CASE("Test create_shares_with_points", "[pss][shares_with_points]") {
    print_separator();
    std::cout << "==== Test create_shares_with_points ====" << std::endl;
    
    ZZ gr_prime = conv<ZZ>(2);
    long positive_integer = 4;
    long gr_degree = 4;

    gr galoisring(gr_prime, positive_integer, gr_degree);

    int members = 5, packed_number = 2, poly_degree = 3;

    packed_shamir::scheme Scheme(members, packed_number, poly_degree, galoisring);

    auto is_invertible = [&](const ZZ_pE& element) {
        try {
            ZZ_pE inv = Inv(element, gr_degree);
            return !IsZero(inv);
        } catch (...) {
            return false;
        }
    };

    vector<ZZ_pE> x_points, y_points;
    for (int i = 0; i < poly_degree + 1; i++) {
        ZZ_pE x, y;
        do {
            x = random_ZZ_pE();
            y = random_ZZ_pE();
        } while (!is_invertible(x) || !is_invertible(y));
        x_points.push_back(x);
        y_points.push_back(y);
    }

    vec_ZZ_pE shares;
    bool success = false;
    int max_attempts = 10;

    for (int attempt = 0; attempt < max_attempts && !success; ++attempt) {
        try {
            shares = Scheme.create_shares_with_points(x_points, y_points);
            success = true;
        } catch (const std::exception& e) {
            if (attempt == max_attempts - 1) {
                FAIL("Failed to create shares after " + std::to_string(max_attempts) + " attempts");
            }
            x_points.clear();
            y_points.clear();
            for (int i = 0; i < poly_degree + 1; i++) {
                ZZ_pE x, y;
                do {
                    x = random_ZZ_pE();
                    y = random_ZZ_pE();
                } while (!is_invertible(x) || !is_invertible(y));
                x_points.push_back(x);
                y_points.push_back(y);
            }
        }
    }

    REQUIRE(success);
    REQUIRE(shares.length() == members);

    ZZ_pEX interpolated_poly;
    vec_ZZ_pE x_vec, y_vec;
    x_vec.SetLength(members);
    y_vec.SetLength(members);

    for (int i = 0; i < members; i++) {
        x_vec[i] = Scheme.alpha_set[i];
        y_vec[i] = shares[i];
    }

    interpolate_for_GR(interpolated_poly, x_vec, y_vec, galoisring.p, galoisring.k, galoisring.r);

    for (int i = 0; i < members; i++) {
        ZZ_pE y = eval(interpolated_poly, Scheme.alpha_set[i]);
        REQUIRE(y == shares[i]);
    }

    std::cout << "Test completed successfully." << std::endl;
}