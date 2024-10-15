#ifndef MULT_GATE_H
#define MULT_GATE_H

#include <vector>
#include <assert.h>

#include "gate.h"

namespace dp {
  class MultGate : public Gate {
  public:
    MultGate() {}; // for the padding gates below
    MultGate(std::shared_ptr<Gate> left, std::shared_ptr<Gate> right) {
      mLeft = left;
      mRight = right;
    };

    FF GetMu() {
      if ( !mLearned ) {
	      throw std::invalid_argument("Run multiplication protocol first");
      }
      return mMu;
    };

    void SetDummyLambda(FF lambda) {
      mLambda = lambda;
      mLambdaSet = true;
    };

    FF GetDummyLambda() {
      if ( !mLambdaSet ) {
	      throw std::invalid_argument("Lambda is not set in this multiplication gate");
      }
      return mLambda;
    };

    void SetIndvShrLambda(Shr indv_shr) {
      mIndvShrLambdaC = indv_shr;
      mIndvShrLambdaCSet = true;
    }

    Shr GetIndvShrLambda() {
      if ( !mIndvShrLambdaCSet )
	      throw std::invalid_argument("IndvShrLambda is not set in this multiplication gate");
      return mIndvShrLambdaC;
    }

    //FF GetDn07Share() {
    //  if ( !mDn07Set )
	  //    throw std::invalid_argument("Dn07 shares is not set in this multiplication gate");
    //  return mDn07Share;
    //}

    FF GetClear() {
      if ( !mEvaluated ) {
	      mClear = mLeft->GetClear() * mRight->GetClear();
	      mEvaluated = true;
      }
      return mClear;
    }

  private:
  };

  // Used for padding batched multiplications
  class PadMultGate : public MultGate {
  public:
    PadMultGate() : MultGate() { mIsPadding = true; };

    FF GetMu() override { return FF(0); }
    FF GetDummyLambda() override { return FF(0); }

    // Just a technicality needed to make the parents of this gate be
    // itself when instantiated
    void UpdateParents(std::shared_ptr<Gate> left, std::shared_ptr<Gate> right) {
      mLeft = left;
      mRight = right;
    }

  private:
  };
    
  class MultBatch {
  public:
    MultBatch(std::size_t batch_size_l, std::size_t batch_size_m) : mBatchSize_l(batch_size_l), mBatchSize_m(batch_size_m) {
      mBatchSize = mBatchSize_l*mBatchSize_m;
      mMultGatesPtrs.reserve(mBatchSize);
      Scheme_m1 = packed_shamir::scheme(mParties, mBatchSize_m, mBatchSize_m-1,gring);
    };

    // Adds a new mult_gate to the batch. It cannot add more gates than
    // the batch_size
    void Append(std::shared_ptr<MultGate> mult_gate) {
      if ( mMultGatesPtrs.size() == mBatchSize )
	      throw std::invalid_argument("Trying to batch more than batch_size gates");
      mMultGatesPtrs.emplace_back(mult_gate); };

    // For testing purposes: sets the required preprocessing for this
    // batch to be just constant shares
    void _DummyPrep(FF lambda_A, FF lambda_B, FF lambda_C) {
      if ( mMultGatesPtrs.size() != mBatchSize )
	      throw std::invalid_argument("The number of mult gates does not match the batch size");
      Shr lambda_A_ = conv<Shr>(lambda_A);
      Shr lambda_B_ = conv<Shr>(lambda_B);
      Shr lambda_C_ = conv<Shr>(lambda_C);

      mPackedShrLambdaA = lambda_A_;
      mPackedShrLambdaB = lambda_B_;
      mPackedShrDeltaC = lambda_A_ * lambda_B_ - lambda_C_;
    };

    void _DummyPrep() {
      _DummyPrep(FF(0), FF(0), FF(0));
    };
    
    //TODO: dont use beaver triple, implement semi_honest version first
    // Generates the preprocessing from the lambdas of the inputs
    void PrepFromDummyLambdas() {
      vector<FF> lambda_A;
      vector<FF> lambda_B;
      vector<FF> delta_C;

      for (std::size_t i = 0; i < mBatchSize; i++) {
	      auto l_A = mMultGatesPtrs[i]->GetLeft()->GetDummyLambda();
	      auto l_B = mMultGatesPtrs[i]->GetRight()->GetDummyLambda();
	      auto l_C = mMultGatesPtrs[i]->GetDummyLambda();
	      auto d_C = l_A * l_B - l_C;

	      lambda_A.emplace_back(l_A);
	      lambda_B.emplace_back(l_B);
	      delta_C.emplace_back(d_C);
      }

      vec_ZZ_pE lambda_A_pE;
      vec_ZZ_pE lambda_B_pE;
      vec_ZZ_pE delta_C_pE;
      lambda_A_pE.SetLength(mBatchSize_m);
      lambda_B_pE.SetLength(mBatchSize_m);
      delta_C_pE.SetLength(mBatchSize_m);
      for(long i=0; i<mBatchSize_m; i++){
        vector<long> temp1;
        vector<long> temp2;
        vector<long> temp3;
        for(long j=0; j<mBatchSize_l; j++){
          temp1[j] = conv<long>(lambda_A[j+i*mBatchSize_l]);
          temp2[j] = conv<long>(lambda_B[j+i*mBatchSize_l]);
          temp3[j] = conv<long>(delta_C[j+i*mBatchSize_l]);
        }
        rmfe.set_input(temp1);
        rmfe.RMFE_GR_PHI();
        vector<long> result1 = rmfe.get_result();
        lambda_A_pE[i] = long2ZZpE(result1);

        rmfe.set_input(temp2);
        rmfe.RMFE_GR_PHI();
        vector<long> result2 = rmfe.get_result();
        lambda_B_pE[i] = long2ZZpE(result2);

        rmfe.set_input(temp3);
        rmfe.RMFE_GR_PHI();
        vector<long> result3 = rmfe.get_result();
        delta_C_pE[i] = long2ZZpE(result3);
      }
      
      // Using deg = BatchSize-1 ensures there's no randomness involved
      //auto poly_A = scl::details::EvPolyFromSecretsAndDegree(lambda_A, mBatchSize-1, mPRG);
      //mPackedShrLambdaA = poly_A.Evaluate(FF(mID));
      vec_ZZ_pE shares_A = Scheme_m1.create_shares(lambda_A_pE);
      mPackedShrLambdaA = shares_A[mID];
	
      //auto poly_B = scl::details::EvPolyFromSecretsAndDegree(lambda_B, mBatchSize-1, mPRG);
      //mPackedShrLambdaB = poly_B.Evaluate(FF(mID));
      vec_ZZ_pE shares_B = Scheme_m1.create_shares(lambda_B_pE);
      mPackedShrLambdaA = shares_B[mID];

      //auto poly_C = scl::details::EvPolyFromSecretsAndDegree(delta_C, mBatchSize-1, mPRG);
      //mPackedShrDeltaC = poly_C.Evaluate(FF(mID));
      vec_ZZ_pE shares_C = Scheme_m1.create_shares(delta_C_pE);
      mPackedShrDeltaC = shares_C[mID];
    }

    // For cleartext evaluation: calls GetClear on all its gates to
    // populate their mClear. This could return a vector with these
    // values but we're not needing them
    void GetClear() {
      for (auto gate : mMultGatesPtrs) { gate->GetClear(); }
    }

    // Determines whether the batch is full
    bool HasRoom() { return mMultGatesPtrs.size() < mBatchSize; }
    
    // For fetching mult gates
    std::shared_ptr<MultGate> GetMultGate(std::size_t idx) { return mMultGatesPtrs[idx]; }

    // Set network parameters for evaluating the protocol. This is not
    // part of the creation of the batch since sometimes we just want
    // to evaluate in the clear and this won't be needed
    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
    }

    // First step of the protocol where P1 sends the packed shares of
    // the mu of the inputs
    void P1Sends();

    // The parties receive the packed shares of the mu's and store
    // them
    void PartiesReceive();

    // The parties compute locally the necessary Beaver linear
    // combination and send the resulting share back to P1
    void PartiesSend();

    // P1 receives the shares from the parties, reconstructs the mu
    // of the outputs in the batch, and updates these accordingly
    void P1Receives();

    void RunProtocol() {
      P1Sends();
      PartiesReceive();
      PartiesSend();
      P1Receives();
    }

    void SetPreprocessing(Shr shr_lambda_A, Shr shr_lambda_B, Shr shr_delta_C) {
      mPackedShrLambdaA = shr_lambda_A;
      mPackedShrLambdaB = shr_lambda_B;
      mPackedShrDeltaC = shr_delta_C;
    }


  private:
    std::size_t mBatchSize;
    std::size_t mBatchSize_l;
    std::size_t mBatchSize_m;

    packed_shamir::scheme Scheme_m1;

    // The mult gates that are part of this batch
    vec<std::shared_ptr<MultGate>> mMultGatesPtrs;

    // The packed sharings associated to this batch
    Shr mPackedShrLambdaA;
    Shr mPackedShrLambdaB;
    Shr mPackedShrDeltaC;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    // Intermediate-protocol
    scl::PRG mPRG;
    Shr mPackedShrMuA;    // Shares of mu_alpha
    Shr mPackedShrMuB;    // Shares of mu_beta
  };

  // Basically a collection of batches
  class MultLayer {
  public:
    MultLayer(std::size_t batch_size_l, std::size_t batch_size_m) : mBatchSize_l(batch_size_l),  mBatchSize_m(batch_size_m){
      mBatchSize = mBatchSize_l*mBatchSize_m;
      auto first_batch = std::make_shared<MultBatch>(mBatchSize_l, mBatchSize_m);
      // Append a first batch
      mBatches.emplace_back(first_batch);
    };

    // Adds a new mult_gate to the layer. It checks if the current
    // batch is full and if so creates a new one.
    void Append(std::shared_ptr<MultGate> mult_gate) {
      auto current_batch = mBatches.back(); // accessing last elt
      if ( current_batch->HasRoom() ) {
	      current_batch->Append(mult_gate);
      } else {
	      auto new_batch = std::make_shared<MultBatch>(mBatchSize_l, mBatchSize_m);
	      new_batch->Append(mult_gate);
	      mBatches.emplace_back(new_batch);
      }
    }

    std::shared_ptr<MultBatch> GetMultBatch(std::size_t idx) { return mBatches[idx]; }

    // Pads the current batch if necessary
    void Close() {
      auto padding_gate = std::make_shared<PadMultGate>();
      padding_gate->UpdateParents(padding_gate, padding_gate);

      assert(padding_gate->GetMu() == FF(0));
      assert(padding_gate->GetLeft()->GetMu() == FF(0));
      assert(padding_gate->GetRight()->GetMu() == FF(0));      

      auto last_batch = mBatches.back(); // accessing last elt
      while ( last_batch->HasRoom() ) {
	      last_batch->Append(padding_gate);
      } 
      // assert(last_batch->HasRoom() == false); // PASSES
    }

    // For testing purposes: sets the required preprocessing for each
    // batch to be just 0 shares
    void _DummyPrep(FF lambda_A, FF lambda_B, FF lambda_C) {
      for (auto batch : mBatches) batch->_DummyPrep(lambda_A, lambda_B, lambda_C);
    }
    void _DummyPrep() {
      for (auto batch : mBatches) batch->_DummyPrep();
    }
    void PrepFromDummyLambdas() {
      for (auto batch : mBatches) batch->PrepFromDummyLambdas();
    }

    void ClearEvaluation() {
      for (auto batch : mBatches) batch->GetClear();
    }

    // Set network parameters for evaluating the protocol. This is not
    // part of the creation of the layer since sometimes we just want
    // to evaluate in the clear and this won't be needed
    // To be used after the layer is closed
    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
      for (auto batch : mBatches) batch->SetNetwork(network, id);
    }

    // For testing purposes, to avoid blocking
    void P1Sends() { for (auto batch : mBatches) batch->P1Sends(); }
    void PartiesReceive() { for (auto batch : mBatches) batch->PartiesReceive(); }
    void PartiesSend() { for (auto batch : mBatches) batch->PartiesSend(); }
    void P1Receives() { for (auto batch : mBatches) batch->P1Receives(); }
    
    // Metrics
    std::size_t GetSize() { return mBatches.size(); }

  private:
    vec<std::shared_ptr<MultBatch>> mBatches;
    std::size_t mBatchSize;
    std::size_t mBatchSize_l;
    std::size_t mBatchSize_m;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    friend class Circuit;
  };
} // namespace tp

#endif  // MULT_GATE_H
