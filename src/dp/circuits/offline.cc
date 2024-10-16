#include "dp/circuits.h"

namespace dp {
  void Circuit::_DummyPrep(FF lambda) {
    for (auto input_gate : mInputGates) {
      input_gate->SetLambda(lambda);
    }
    for (std::size_t layer = 0; layer < GetDepth(); layer++) {
      for (auto mult_gate : mFlatMultLayers[layer]) {
	      mult_gate->SetDummyLambda(lambda);
      }
    }
    for (auto output_gate : mOutputGates) {
      output_gate->GetDummyLambda(); //populate outputs and add wires
    }

    // Now update packed sharings
    for (auto input_layer : mInputLayers) {
      for (auto input_batch : input_layer.mBatches) {
        vec<FF> lambda_A_FF;
	      Vec lambda_A;
	      lambda_A.reserve(mBatch_m);
        lambda_A_FF.reserve(mBatchSize);
	      for (std::size_t i = 0; i < mBatchSize; i++) {
	        lambda_A_FF.emplace_back(input_batch->GetInputGate(i)->GetDummyLambda());
	      }

        for(std::size_t i=0; i<mBatch_m; i++){
          vector<long> temp;
          for(std::size_t j=0; j<mBatch_l; j++){
            temp[j] = conv<long>(lambda_A_FF[j+i*mBatch_l]);
          }
          rmfe.set_input(temp);
          rmfe.RMFE_GR_PHI();
          vector<long> result = rmfe.get_result();
          lambda_A[i] = long2ZZpE(result);
        }
	      //scl::PRG prg;
	      //auto poly = scl::details::EvPolyFromSecretsAndDegree(lambda_A, mBatchSize-1, prg);
	      //Vec new_shares = scl::details::SharesFromEvPoly(poly, mParties);

        vec_ZZ_pE new_shares = Scheme_m1.create_shares(lambda_A);

	      input_batch->SetPreprocessing(new_shares[mID]);
      }
    }
    for (auto output_layer : mOutputLayers) {
      for (auto output_batch : output_layer.mBatches) {
	      //Vec lambda_A;
	      //lambda_A.reserve(mBatch_m);
	      //for (std::size_t i = 0; i < mBatch_m; i++) {
	      //  lambda_A.emplace_back(output_batch->GetOutputGate(i)->GetDummyLambda());
	      //}
	      //scl::PRG prg;
	      //auto poly = scl::details::EvPolyFromSecretsAndDegree(lambda_A, mBatchSize-1, prg);
	      //Vec new_shares = scl::details::SharesFromEvPoly(poly, mParties);
        //vec_ZZ_pE new_shares = Scheme_m1.create_shares(lambda_A);

        vec<FF> lambda_A_FF;
	      Vec lambda_A;
	      lambda_A.reserve(mBatch_m);
        lambda_A_FF.reserve(mBatchSize);
	      for (std::size_t i = 0; i < mBatchSize; i++) {
	        lambda_A_FF.emplace_back(output_batch->GetOutputGate(i)->GetDummyLambda());
	      }

        for(std::size_t i=0; i<mBatch_m; i++){
          vector<long> temp;
          for(std::size_t j=0; j<mBatch_l; j++){
            temp[j] = conv<long>(lambda_A_FF[j+i*mBatch_l]);
          }
          rmfe.set_input(temp);
          rmfe.RMFE_GR_PHI();
          vector<long> result = rmfe.get_result();
          lambda_A[i] = long2ZZpE(result);
        }

        vec_ZZ_pE new_shares = Scheme_m1.create_shares(lambda_A);

	      output_batch->SetPreprocessing(new_shares[mID]);
      }
    }
    // Mults
    for (auto mult_layer : mMultLayers) {
      for (auto mult_batch : mult_layer.mBatches) {
	      vec<FF> lambda_A;
	      vec<FF> lambda_B;
	      vec<FF> lambda_C;
	      lambda_A.reserve(mBatchSize);
	      lambda_B.reserve(mBatchSize);
	      lambda_C.reserve(mBatchSize);
	      for (std::size_t i = 0; i < mBatchSize; i++) {
	        lambda_A.emplace_back(mult_batch->GetMultGate(i)->GetLeft()->GetDummyLambda());
	        lambda_B.emplace_back(mult_batch->GetMultGate(i)->GetRight()->GetDummyLambda());
	        lambda_C.emplace_back(mult_batch->GetMultGate(i)->GetDummyLambda());
	      }

        Vec lambda_A_;
        Vec lambda_B_;
        Vec lambda_C_;
	      lambda_A_.reserve(mBatch_m);
        lambda_B_.reserve(mBatch_m);
        lambda_C_.reserve(mBatch_m);

        for(std::size_t i=0; i<mBatch_m; i++){
          vector<long> temp1;
          vector<long> temp2;
          vector<long> temp3;
          for(std::size_t j=0; j<mBatch_l; j++){
            temp1[j] = conv<long>(lambda_A[j+i*mBatch_l]);
            temp2[j] = conv<long>(lambda_B[j+i*mBatch_l]);
            temp3[j] = conv<long>(lambda_C[j+i*mBatch_l]);
          }
          rmfe.set_input(temp1);
          rmfe.RMFE_GR_PHI();
          vector<long> result1 = rmfe.get_result();
          lambda_A_[i] = long2ZZpE(result1);
          rmfe.set_input(temp2);
          rmfe.RMFE_GR_PHI();
          vector<long> result2 = rmfe.get_result();
          lambda_B_[i] = long2ZZpE(result2);
          rmfe.set_input(temp3);
          rmfe.RMFE_GR_PHI();
          vector<long> result3 = rmfe.get_result();
          lambda_C_[i] = long2ZZpE(result3);
        }

	      //scl::PRG prg;
	      //auto poly_A = scl::details::EvPolyFromSecretsAndDegree(lambda_A, mBatchSize-1, prg);
	      //auto poly_B = scl::details::EvPolyFromSecretsAndDegree(lambda_B, mBatchSize-1, prg);
	      //auto poly_C = scl::details::EvPolyFromSecretsAndDegree(lambda_C, mBatchSize-1, prg);
	      vec_ZZ_pE new_shares_A = Scheme_m1.create_shares(lambda_A_);
	      vec_ZZ_pE new_shares_B = Scheme_m1.create_shares(lambda_B_);
	      vec_ZZ_pE new_shares_C = Scheme_m1.create_shares(lambda_C_);

	      mult_batch->SetPreprocessing(new_shares_A[mID], new_shares_B[mID], \
				     new_shares_A[mID] * new_shares_B[mID] - new_shares_C[mID]);
      }
    }
      
      
  }

    // Populates each batch with dummy preprocessing (all zeros)
  void Circuit::_DummyPrep() {
      for (auto input_layer : mInputLayers) input_layer._DummyPrep();
      for (auto mult_layer : mMultLayers) mult_layer._DummyPrep();
      for (auto output_layer : mOutputLayers) output_layer._DummyPrep();
    }

    // Set all settable lambdas to a constant
    void Circuit::SetDummyLambdas(FF lambda) {
      // Output wires of multiplications
      for (std::size_t layer = 0; layer < GetDepth(); layer++) {
	      for (auto mult_gate : mFlatMultLayers[layer]) {
	        mult_gate->SetDummyLambda(lambda);
	      }
      }
      // Output wires of input gates
      for (auto input_gate : mInputGates) {
	      input_gate->SetLambda(lambda);
      }
    }

    // Populates the lambdas of addition and output gates
    void Circuit::PopulateDummyLambdas() {
      // Addition gates
      for (auto add_gate : mAddGates) {
	      (void)add_gate->GetDummyLambda();
      }
      // Output gates
      for (auto output_gate : mOutputGates) {
	      (void)output_gate->GetDummyLambda();
      }      
    }

    // Sets the lambdas for the output wires of addition and output
    // gates based on the lambdas for multiplications and input gates
    void Circuit::PrepFromDummyLambdas() {
      for (auto input_layer : mInputLayers) input_layer.PrepFromDummyLambdas();
      for (auto output_layer : mOutputLayers) output_layer.PrepFromDummyLambdas();
      for (auto mult_layer : mMultLayers) mult_layer.PrepFromDummyLambdas();
    }

    void Circuit::GenCorrelator() {
      if ( !mIsNetworkSet )
	      throw std::invalid_argument("Cannot set correlator without setting a network first");
      
      std::size_t n_ind_shares = GetNInputs() + GetSize();
      std::size_t n_mult_batches = GetNMultBatches();
      std::size_t n_inout_batches = GetNInputBatches() + GetNOutputBatches();

      mCorrelator = Correlator(n_ind_shares, n_mult_batches, n_inout_batches, mBatchSize, mBatch_l, mBatch_m, Scheme_nm,Scheme_n1,Scheme_t,Scheme_m1);
      mCorrelator.SetNetwork(mNetwork, mID);
      mCorrelator.PrecomputeEij();
    }

    // Populates the mappings in the correlator so that the
    // F.I. preprocessing is mapped to different gates/batches
    void Circuit::MapCorrToCircuit() {
      // Individual shares
      for (auto input_gate : mInputGates) {
	      mCorrelator.PopulateIndvShrs(input_gate);
      }
      for (std::size_t layer = 0; layer < GetDepth(); layer++) {
	      for (auto mult_gate : mFlatMultLayers[layer]) {
	        mCorrelator.PopulateIndvShrs(mult_gate);
	      }
      }
      for (auto add_gate : mAddGates) {
	      mCorrelator.PopulateIndvShrs(add_gate);
      }
      for (auto output_gate : mOutputGates) {
	      mCorrelator.PopulateIndvShrs(output_gate);
      }

      // Input batches
      for (auto input_layer : mInputLayers) {
	      for (auto input_batch : input_layer.mBatches) {
	        mCorrelator.PopulateInputBatches(input_batch);
	      }
      }
      // Output batches
      for (auto output_layer : mOutputLayers) {
	      for (auto output_batch : output_layer.mBatches) {
	        mCorrelator.PopulateOutputBatches(output_batch);
	      }
      }
      // Mult batches
      for (auto mult_layer : mMultLayers) {
	      for (auto mult_batch : mult_layer.mBatches) {
	        mCorrelator.PopulateMultBatches(mult_batch);
	      }
      }
    }

    // Prep inputs & outputs
    void Circuit::PrepMultPartiesSendP1() {
      for (auto mult_layer : mMultLayers) {
	      for (auto mult_batch : mult_layer.mBatches) {
	        mCorrelator.PrepMultPartiesSendP1(mult_batch);
	      }
      }
    }
    void Circuit::PrepMultP1ReceivesAndSends() {
      for (auto mult_layer : mMultLayers) {
	      for (auto mult_batch : mult_layer.mBatches) {
	        (void) mult_batch;
	        mCorrelator.PrepMultP1ReceivesAndSends();
	      } 
      }      
    }
    void Circuit::PrepMultPartiesReceive() {
      for (auto mult_layer : mMultLayers) {
	      for (auto mult_batch : mult_layer.mBatches) {
	        mCorrelator.PrepMultPartiesReceive(mult_batch);
	      }
      }
    }

    void Circuit::PrepIOPartiesSendOwner() {
      for (auto input_layer : mInputLayers) {
	      for (auto input_batch : input_layer.mBatches) {
	        mCorrelator.PrepInputPartiesSendOwner(input_batch);
	      }
      }
      for (auto output_layer : mOutputLayers) {
	      for (auto output_batch : output_layer.mBatches) {
	        mCorrelator.PrepOutputPartiesSendOwner(output_batch);
	      }
      }
    }

    void Circuit::PrepIOOwnerReceives() {
      for (auto input_layer : mInputLayers) {
	      for (auto input_batch : input_layer.mBatches) {
	        mCorrelator.PrepInputOwnerReceives(input_batch);
	      } 
      }
      for (auto output_layer : mOutputLayers) {
	      for (auto output_batch : output_layer.mBatches) {
	        mCorrelator.PrepOutputOwnerReceives(output_batch);
	      }
      }
    }
  
} // namespace tp
