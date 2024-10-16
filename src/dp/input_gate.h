#ifndef INP_GATE_H
#define INP_GATE_H

#include <vector>
#include <assert.h>

#include "gate.h"

namespace dp {
  class InputGate : public Gate {
  public:
    InputGate(std::size_t owner_id) : mOwnerID(owner_id) {
      // TODO sample at random, interactively. Get from correlator
      mIndvShrLambdaC = Shr(0);
    }

    // Set cleartext inputs, for the case of cleartext evaluation
    void ClearInput(FF input) {
      mClear = input;
      mEvaluated = true;
    }

    std::size_t GetOwner() {
      return mOwnerID;
    }

    FF GetMu() {
      if ( !mLearned )
	      throw std::invalid_argument("P1 hasn't learned this value yet");
      return mMu;
    }

    void SetLambda(FF lambda) {
      mLambda = lambda;
      mLambdaSet = true;
    }

    FF GetDummyLambda() {
      if ( !mLambdaSet )
	      throw std::invalid_argument("Lambda is not set in this input gate");
      return mLambda;
    }

    void SetIndvShrLambda(Shr indv_shr) {
      mIndvShrLambdaC = indv_shr;
      mIndvShrLambdaCSet = true;
    }

    Shr GetIndvShrLambda() {
      if ( !mIndvShrLambdaCSet )
	      throw std::invalid_argument("IndvShrLambda is not set in this input gate");
      return mIndvShrLambdaC;
    }

    //FF GetDn07Share() {
    //  if ( !mDn07Set )
	  //    throw std::invalid_argument("Dn07 shares is not set in this input gate");
    //  return mDn07Share;
    //}
    
    FF GetClear() {
      if ( !mEvaluated )
	      throw std::invalid_argument("This input has not been provided yet");
      return mClear;
    };

    // Protocol-related
    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
    }

    void _DummyPrep(FF lambda) {
      if (mID == mOwnerID) mLambda = lambda;
    }
    void _DummyPrep() {
      _DummyPrep(FF(0));
    }

    void SetInput(FF input) {
      if ( mID == mOwnerID ) mValue = input;
    }

    void OwnerSendsP1() {
      if (mID == mOwnerID) mNetwork->Party(0)->Send(mValue - mLambda);
    }

    void P1Receives() {
      if (mID == 0) {
	      mNetwork->Party(mOwnerID)->Recv(mMu);
	      mLearned = true;
      }
    }

    FF GetValue() { return mValue; }

  private:
    // ID of the party who owns this gate
    std::size_t mOwnerID;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    // Protocol-specific
    FF mLambda; // Lambda, learned by owner
    FF mValue; // Actual input, known by owner
  };

  // Used for padding batched inputs
  class PadInputGate : public InputGate {
  public:
    PadInputGate(std::size_t owner_id) : InputGate(owner_id) { mIsPadding = true; }

    FF GetMu() override { return FF(0); }
    FF GetDummyLambda() override { return FF(0); }
    Shr GetIndvShrLambda() override { return Shr(0); }

  private:
  };

  class InputBatch {
  public:
    InputBatch(std::size_t owner_id, std::size_t batch_size_l, std::size_t batch_size_m) : mOwnerID(owner_id), mBatchSize_l(batch_size_l),  mBatchSize_m(batch_size_m){
      mBatchSize = mBatchSize_l*mBatchSize_m;
      mInputGatesPtrs.reserve(mBatchSize); //TODO:should add one layer
      Scheme_m1 = packed_shamir::scheme(mParties, mBatchSize_m, mBatchSize_m-1,gring);
    };

    // Adds a new input_gate to the batch. It cannot add more gates than
    // the batch_size
    void Append(std::shared_ptr<InputGate> input_gate) {
      if ( mInputGatesPtrs.size() == mBatchSize )
	      throw std::invalid_argument("Trying to batch more than batch_size gates");
      if ( input_gate->GetOwner() != mOwnerID )
	      throw std::invalid_argument("Owner IDs do not match");
      mInputGatesPtrs.emplace_back(input_gate); }

    // For testing purposes: sets the required preprocessing for this
    // batch to be constant shares
    void _DummyPrep(FF lambda) {
      if ( mInputGatesPtrs.size() != mBatchSize )
	      throw std::invalid_argument("The number of input gates does not match the batch size");
      Shr lambda_ = conv<Shr>(lambda);
      mPackedShrLambda = lambda_;
      for (auto input_gate : mInputGatesPtrs) input_gate->_DummyPrep(lambda);
    }
    void _DummyPrep() {
      _DummyPrep(FF(0));
    }


    std::size_t GetOwner() {
      return mOwnerID;
    };

                                            
    // Generates the preprocessing from the lambdas of the inputs
    void PrepFromDummyLambdas() {
      vector<FF> lambda;//TODO:dont understand. exist problems in turbopack i think

      for (std::size_t i = 0; i < mBatchSize; i++) {
	      lambda.emplace_back(mInputGatesPtrs[i]->GetDummyLambda());
      }

      vec_ZZ_pE lambda_;
      lambda_.SetLength(mBatchSize_m);
      for(long i=0; i<mBatchSize_m; i++){
        vector<long> temp;
        for(long j=0; j<mBatchSize_l; j++){
          temp[j] = conv<long>(lambda[j+i*mBatchSize_l]);
        }
        rmfe.set_input(temp);
        rmfe.RMFE_GR_PHI();
        vector<long> result = rmfe.get_result();
        lambda_[i] = long2ZZpE(result);
      }
      // Using deg = BatchSize-1 ensures there's no randomness involved
      //auto poly = scl::details::EvPolyFromSecretsAndDegree(lambda, mBatchSize-1, mPRG);
      //Vec shares = scl::details::SharesFromEvPoly(poly, mParties);
      vec_ZZ_pE shares = Scheme_m1.create_shares(lambda_);

      mPackedShrLambda = shares[mID];
    }


    // For cleartext evaluation: calls GetClear on all its gates to
    // populate their mClear. This could return a vector with these
    // values but we're not needing them
    void GetClear() {
      for (auto gate : mInputGatesPtrs) { gate->GetClear(); }
    }

    // Determines whether the batch is full
    bool HasRoom() { return mInputGatesPtrs.size() < mBatchSize; }
    
    // For fetching input gates
    std::shared_ptr<InputGate> GetInputGate(std::size_t idx) { return mInputGatesPtrs[idx]; }

    // Set network parameters for evaluating the protocol. This is not
    // part of the creation of the batch since sometimes we just want
    // to evaluate in the clear and this won't be needed
    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
      for (auto input_gate : mInputGatesPtrs) input_gate->SetNetwork(network, id);
    }

    void SetPreprocessing(Shr packed_shr_lambda) { mPackedShrLambda = packed_shr_lambda; }

    Shr GetPackedShrLambda() { return mPackedShrLambda; }


  private:
    // ID of the party who owns this batch
    std::size_t mOwnerID;

    std::size_t mBatchSize; // mBatchSize = ml
    std::size_t mBatchSize_m;
    std::size_t mBatchSize_l;


    packed_shamir::scheme Scheme_m1; 
    // The input gates that are part of this batch
    vec<std::shared_ptr<InputGate>> mInputGatesPtrs;

    // The packed sharings associated to this batch
    Shr mPackedShrLambda;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    scl::PRG mPRG;
  };

  // Basically a collection of batches
  class InputLayer {
  public:
    InputLayer(std::size_t owner_id, std::size_t batch_size_l, std::size_t batch_size_m) : mOwnerID(owner_id), mBatchSize_l(batch_size_l),  mBatchSize_m(batch_size_m) {
      mBatchSize = mBatchSize_l*mBatchSize_m;
      auto first_batch = std::make_shared<InputBatch>(mOwnerID, mBatchSize_l, mBatchSize_m);
      // Append a first batch
      mBatches.emplace_back(first_batch);
    };

    // Adds a new input_gate to the layer. It checks if the current
    // batch is full and if so creates a new one.
    void Append(std::shared_ptr<InputGate> input_gate) {
      auto current_batch = mBatches.back(); // accessing last elt
      if ( current_batch->HasRoom() ) {
	      current_batch->Append(input_gate);
      } else {
	      auto new_batch = std::make_shared<InputBatch>(mOwnerID, mBatchSize_l, mBatchSize_m);
	      new_batch->Append(input_gate);
	      mBatches.emplace_back(new_batch);
      }
    }

    std::shared_ptr<InputBatch> GetInputBatch(std::size_t idx) { return mBatches[idx]; }

    // Pads the current batch if necessary
    void Close() {
      auto padding_gate = std::make_shared<PadInputGate>(mOwnerID);
      assert(padding_gate->GetMu() == FF(0));

      auto last_batch = mBatches.back(); // accessing last elt
      while ( last_batch->HasRoom() ) {
	      last_batch->Append(padding_gate);
      } 
      // assert(last_batch->HasRoom() == false); // PASSES
    }

    // For testing purposes: sets the required preprocessing for each
    // batch to be just 0 shares
    void _DummyPrep(FF lambda) {
      for (auto batch : mBatches) batch->_DummyPrep(lambda);
    }
    void _DummyPrep() {
      for (auto batch : mBatches) batch->_DummyPrep();
    }
    void PrepFromDummyLambdas() {
      for (auto batch : mBatches) batch->PrepFromDummyLambdas();
    }

    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
      for (auto batch : mBatches) batch->SetNetwork(network, id);
    }

    void ClearEvaluation() {
      for (auto batch : mBatches) batch->GetClear();
    }

    // Metrics
    std::size_t GetSize() { return mBatches.size(); }

  private:
    std::size_t mOwnerID;
    vec<std::shared_ptr<InputBatch>> mBatches;
    std::size_t mBatchSize;
    std::size_t mBatchSize_m;
    std::size_t mBatchSize_l;

    // Network-related
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    friend class Circuit;
  };

} // namespace tp

#endif  // INP_GATE_H
