#include <catch2/catch.hpp>
#include <iostream>

#include "dp/input_gate.h"
#include "dp/mult_gate.h"

TEST_CASE("Mult") {
  std::size_t threshold = 4; // has to be even
  std::size_t batch_size_m = (threshold + 2)/2;
  std::size_t batch_size_l = 2;
  std::size_t n_parties = threshold + 2*(batch_size_m - 1) + 1;
  auto networks = scl::Network::CreateFullInMemory(n_parties);

  dp::FF lambda_x(239);
  dp::FF lambda_y(-3421);
  dp::FF lambda_z(942582);

  std::vector<std::shared_ptr<dp::InputGate>> x_gates;
  std::vector<std::shared_ptr<dp::InputGate>> y_gates;
  x_gates.reserve(n_parties);
  y_gates.reserve(n_parties);

  for (std::size_t i = 0; i < n_parties; i++) {
    x_gates.emplace_back(std::make_shared<dp::InputGate>(0));
    y_gates.emplace_back(std::make_shared<dp::InputGate>(0));
  }

  dp::FF X(21321);
  dp::FF Y(-3421);

  for (std::size_t i = 0; i < n_parties; i++) {
    x_gates[i]->_SetDummyMu(X - lambda_x);
    y_gates[i]->_SetDummyMu(Y - lambda_y);
  }  
  
  SECTION("MultBatch"){
    std::vector<std::shared_ptr<dp::MultGate>> z_gates;
    z_gates.reserve(n_parties);

    for (std::size_t i = 0; i < n_parties; i++) {
      z_gates.emplace_back(std::make_shared<dp::MultGate>( x_gates[i], y_gates[i] ));
    }

    // Setting batches
    std::vector<std::shared_ptr<dp::MultBatch>> batches;
    batches.reserve(n_parties);

    for (std::size_t i = 0; i < n_parties; i++) {
      batches.emplace_back(std::make_shared<dp::MultBatch>(1,1));
      batches[i]->Append(z_gates[i]);
      batches[i]->SetNetwork(std::make_shared<scl::Network>(networks[i]), i);
      batches[i]->_DummyPrep(lambda_x, lambda_y, lambda_z);
    }

    for (std::size_t i = 0; i < n_parties; i++) {
      batches[i]->P1Sends(); // Only P1 works here
    }
    
    for (std::size_t i = 0; i < n_parties; i++) {
      batches[i]->PartiesReceive(); 
    }

    for (std::size_t i = 0; i < n_parties; i++) {
      batches[i]->PartiesSend(); 
    }

    REQUIRE(!z_gates[0]->IsLearned()); // P1 hasn't learned anything yet

    for (std::size_t i = 0; i < n_parties; i++) {
      batches[i]->P1Receives(); // Only P1 works here
    }

    REQUIRE(z_gates[0]->IsLearned()); // P1 learned mu already
    REQUIRE(!z_gates[1]->IsLearned()); // But not the other parties

    // Check result
    
    dp::FF Z = X*Y;
    dp::FF mu_z = z_gates[0]->GetMu(); 
    REQUIRE( Z == mu_z + lambda_z ); 
  }

  SECTION("MultLayer"){
    std::vector<std::shared_ptr<dp::MultGate>> z_gates;
    z_gates.reserve(n_parties);

    for (std::size_t i = 0; i < n_parties; i++) {
      z_gates.emplace_back(std::make_shared<dp::MultGate>( x_gates[i], y_gates[i] ));
    }

    // Setting layer
    std::vector<dp::MultLayer> layers;
    layers.reserve(n_parties);

    for (std::size_t i = 0; i < n_parties; i++) {
      auto layer = dp::MultLayer(batch_size_l, batch_size_m);
      layer.Append(z_gates[i]);
      layer.Close();
      layer.SetNetwork(std::make_shared<scl::Network>(networks[i]), i);

      layer._DummyPrep(lambda_x, lambda_y, lambda_z); 

      layers.emplace_back(layer);
    }

    // Run the protocol
    for (std::size_t i = 0; i < n_parties; i++) {
      layers[i].P1Sends(); // Only P1 works here
    }
    
    for (std::size_t i = 0; i < n_parties; i++) {
      layers[i].PartiesReceive(); 
    }

    for (std::size_t i = 0; i < n_parties; i++) {
      layers[i].PartiesSend(); 
    }

    REQUIRE(!z_gates[0]->IsLearned()); // P1 hasn't learned anything yet

    for (std::size_t i = 0; i < n_parties; i++) {
      layers[i].P1Receives(); // Only P1 works here
    }

    REQUIRE(z_gates[0]->IsLearned()); // P1 learned mu already
    REQUIRE(!z_gates[1]->IsLearned()); // But not the other parties

    // Check result
    
    dp::FF Z = X*Y;
    dp::FF mu_z = z_gates[0]->GetMu();
    REQUIRE( Z == mu_z + lambda_z ); 
  }

}