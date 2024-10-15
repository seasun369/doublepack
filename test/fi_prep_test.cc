#include <catch2/catch.hpp>
#include "dp/correlator.h"
#include "dp/fi_prep.h"
#include <memory>
#include "../PSS/pss.h"
#include "../RMFE/RMFE.h"
#include "scl/net/network.h"
#include "scl/net/config.h"

using namespace dp;
using namespace packed_shamir;
using namespace scl;

TEST_CASE("FI Prep Basic Functionality", "[fi_prep]") {
    // Setup
    const size_t n_parties = 5;
    const size_t threshold = 2;
    const size_t batch_size = 2;
    const size_t n_mult_batches = 2;
    const size_t n_in_out_batches = 2;
    const size_t n_ind_shrs = 10;

    // Initialize Galois Ring
    ZZ gr_prime = conv<ZZ>(2);
    long positive_integer = 16;
    long gr_degree = 20;
    gr galoisring(gr_prime, positive_integer, gr_degree);

    auto networks = Network::CreateFullInMemory(n_parties);

    SECTION("GenIndShrsPartiesSend and GenIndShrsPartiesReceive") {
        std::vector<std::unique_ptr<Correlator>> correlators;
        for (size_t i = 0; i < n_parties; ++i) {
            std::cout << "Creating Correlator for party " << i << std::endl;
            std::cout << "n_ind_shrs: " << n_ind_shrs << ", n_mult_batches: " << n_mult_batches 
                      << ", n_in_out_batches: " << n_in_out_batches << ", batch_size: " << batch_size << std::endl;
            
            // Adjust the parameters to ensure d-m >= 0
            auto scheme_nm = scheme(n_parties, batch_size, n_parties - batch_size, galoisring);
            auto scheme_n1 = scheme(n_parties, 1, n_parties - 1, galoisring);
            auto scheme_t = scheme(threshold + 1, batch_size, threshold, galoisring);
            auto scheme_m1 = scheme(batch_size, 1, batch_size - 1, galoisring);
            
            std::cout << "scheme_nm: m = " << scheme_nm.get_m() << ", n = " << scheme_nm.get_n() << ", d = " << scheme_nm.get_d() << std::endl;
            std::cout << "scheme_n1: m = " << scheme_n1.get_m() << ", n = " << scheme_n1.get_n() << ", d = " << scheme_n1.get_d() << std::endl;
            std::cout << "scheme_t: m = " << scheme_t.get_m() << ", n = " << scheme_t.get_n() << ", d = " << scheme_t.get_d() << std::endl;
            std::cout << "scheme_m1: m = " << scheme_m1.get_m() << ", n = " << scheme_m1.get_n() << ", d = " << scheme_m1.get_d() << std::endl;
            
            REQUIRE(scheme_nm.get_m() > 0);
            REQUIRE(scheme_nm.get_n() > 0);
            REQUIRE(scheme_nm.get_d() >= scheme_nm.get_m());
            REQUIRE(scheme_n1.get_m() > 0);
            REQUIRE(scheme_n1.get_n() > 0);
            REQUIRE(scheme_n1.get_d() >= scheme_n1.get_m());
            REQUIRE(scheme_t.get_m() > 0);
            REQUIRE(scheme_t.get_n() > 0);
            REQUIRE(scheme_t.get_d() >= scheme_t.get_m());
            REQUIRE(scheme_m1.get_m() > 0);
            REQUIRE(scheme_m1.get_n() > 0);
            REQUIRE(scheme_m1.get_d() >= scheme_m1.get_m());
            
            correlators.push_back(std::make_unique<Correlator>(
                n_ind_shrs, n_mult_batches, n_in_out_batches, batch_size, 
                batch_size, batch_size, scheme_nm, scheme_n1, scheme_t, scheme_m1));
            correlators.back()->SetNetwork(std::make_shared<scl::Network>(networks[i]), i);
            correlators.back()->SetThreshold(threshold);
        }

        // Run GenIndShrsPartiesSend for all parties
        for (size_t i = 0; i < n_parties; ++i) {
            std::cout << "Running GenIndShrsPartiesSend for party " << i << std::endl;
            try {
                REQUIRE(correlators[i] != nullptr);
                correlators[i]->GenIndShrsPartiesSend();
            } catch (const std::exception& e) {
                std::cerr << "Exception caught for party " << i << ": " << e.what() << std::endl;
                FAIL("Exception thrown in GenIndShrsPartiesSend");
            }
        }

        // Run GenIndShrsPartiesReceive for all parties
        for (size_t i = 0; i < n_parties; ++i) {
            std::cout << "Running GenIndShrsPartiesReceive for party " << i << std::endl;
            try {
                REQUIRE(correlators[i] != nullptr);
                correlators[i]->GenIndShrsPartiesReceive();
            } catch (const std::exception& e) {
                std::cerr << "Exception caught for party " << i << ": " << e.what() << std::endl;
                FAIL("Exception thrown in GenIndShrsPartiesReceive");
            }
        }

        // Check that all parties have the correct number of independent shares
        for (size_t i = 0; i < n_parties; ++i) {
            std::cout << "Checking IndShrsSize for party " << i << std::endl;
            REQUIRE(correlators[i]->GetIndShrsSize() == n_ind_shrs);
        }
    }

    SECTION("GenUnpackedShrPartiesSend and GenUnpackedShrPartiesReceive") {
        std::vector<std::unique_ptr<Correlator>> correlators;
        for (size_t i = 0; i < n_parties; ++i) {
            auto scheme_nm = scheme(n_parties, batch_size, n_parties - batch_size, galoisring);
            auto scheme_n1 = scheme(n_parties, 1, n_parties - 1, galoisring);
            auto scheme_t = scheme(threshold + 1, batch_size, threshold, galoisring);
            auto scheme_m1 = scheme(batch_size, 1, batch_size - 1, galoisring);
            
            correlators.push_back(std::make_unique<Correlator>(
                n_ind_shrs, n_mult_batches, n_in_out_batches, batch_size, 
                batch_size, batch_size, scheme_nm, scheme_n1, scheme_t, scheme_m1));
            correlators.back()->SetNetwork(std::make_shared<scl::Network>(networks[i]), i);
            correlators.back()->SetThreshold(threshold);
        }

        // Run GenUnpackedShrPartiesSend for all parties
        for (auto& corr : correlators) {
            corr->GenUnpackedShrPartiesSend();
        }

        // Run GenUnpackedShrPartiesReceive for all parties
        for (auto& corr : correlators) {
            corr->GenUnpackedShrPartiesReceive();
        }

        // Check that all parties have the correct number of unpacked shares
        for (auto& corr : correlators) {
            REQUIRE(corr->GetUnpackedShrsASize() == batch_size);
            REQUIRE(corr->GetUnpackedShrsBSize() == batch_size);
            for (size_t i = 0; i < batch_size; ++i) {
                REQUIRE(corr->GetUnpackedShrsABatchSize(i) == n_mult_batches);
                REQUIRE(corr->GetUnpackedShrsBBatchSize(i) == n_mult_batches);
            }
        }
    }

    SECTION("GenZeroPartiesSend and GenZeroPartiesReceive") {
        std::vector<std::unique_ptr<Correlator>> correlators;
        for (size_t i = 0; i < n_parties; ++i) {
            auto scheme_nm = scheme(n_parties, batch_size, n_parties - batch_size, galoisring);
            auto scheme_n1 = scheme(n_parties, 1, n_parties - 1, galoisring);
            auto scheme_t = scheme(threshold + 1, batch_size, threshold, galoisring);
            auto scheme_m1 = scheme(batch_size, 1, batch_size - 1, galoisring);
            
            correlators.push_back(std::make_unique<Correlator>(
                n_ind_shrs, n_mult_batches, n_in_out_batches, batch_size, 
                batch_size, batch_size, scheme_nm, scheme_n1, scheme_t, scheme_m1));
            correlators.back()->SetNetwork(std::make_shared<scl::Network>(networks[i]), i);
            correlators.back()->SetThreshold(threshold);
        }

        // Run GenZeroPartiesSend for all parties
        for (auto& corr : correlators) {
            corr->GenZeroPartiesSend();
        }

        // Run GenZeroPartiesReceive for all parties
        for (auto& corr : correlators) {
            corr->GenZeroPartiesReceive();
        }

        // Check that all parties have the correct number of zero shares
        for (auto& corr : correlators) {
            REQUIRE(corr->GetMultBatchFIPrepSize() == n_mult_batches);
            REQUIRE(corr->GetIOBatchFIPrepSize() == n_in_out_batches);
        }
    }

    // Add more sections for other functions as needed
}
