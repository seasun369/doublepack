#ifndef CORRELATOR_H
#define CORRELATOR_H

#include <iostream>
#include <map>
#include <assert.h>

#include "dp.h"
#include "mult_gate.h"
#include "input_gate.h"
#include "output_gate.h"
#include "../PSS/pss.h"
using namespace packed_shamir;

namespace dp {
  struct MultBatchFIPrep {
    // Triples
    Shr mShrA;
    Shr mShrB;
    Shr mShrC;

    // Shares of 0
    Shr mShrO1;
    Shr mShrO2;
    Shr mShrO3;

    MultBatchFIPrep(Shr lambda) {
      mShrA = lambda;
      mShrB = lambda;
      mShrC = lambda*lambda;
      mShrO1 = Shr(0);
      mShrO2 = Shr(0);
      mShrO3 = Shr(0);
    }
    MultBatchFIPrep() {
      mShrA = Shr(0);
      mShrB = Shr(0);
      mShrC = Shr(0);
      mShrO1 = Shr(0);
      mShrO2 = Shr(0);
      mShrO3 = Shr(0);
    }
  };

  struct IOBatchFIPrep {
    // Share of zero
    Shr mShrO;

    IOBatchFIPrep() {
      mShrO = Shr(0);
    }
  };  

  class Correlator {
  public:
    Correlator() {};

    Correlator(std::size_t n_ind_shares, std::size_t n_mult_batches, std::size_t n_inout_batches, std::size_t batch_size, std::size_t batch_size_l, std::size_t batch_size_m, scheme scheme_nm, scheme scheme_n1, scheme scheme_t, scheme scheme_m1) :
      mNIndShrs(n_ind_shares), mNMultBatches(n_mult_batches), mNInOutBatches(n_inout_batches), \
      mCTRIndShrs(0), mCTRMultBatches(0), mCTRInOutBatches(0), mBatchSize(batch_size),mBatch_l(batch_size_l),mBatch_m(batch_size_m),Scheme_nm(scheme_nm),Scheme_n1(scheme_n1),Scheme_t(scheme_t), Scheme_m1(scheme_m1) {

      mIndShrs.reserve(mNIndShrs);
      mMultBatchFIPrep.reserve(mNMultBatches);
      mIOBatchFIPrep.reserve(mNInOutBatches);
    }

    // Generates dummy F.I. preprocessing
    void GenIndShrsDummy(Shr lambda) {
      mIndShrs = std::vector<Shr>(mNIndShrs, lambda);
    }
    void GenMultBatchDummy() {
      mMultBatchFIPrep = std::vector<MultBatchFIPrep>(mNMultBatches, MultBatchFIPrep());      
    }
    void GenIOBatchDummy() {
      mIOBatchFIPrep = std::vector<IOBatchFIPrep>(mNInOutBatches, IOBatchFIPrep());
    }

    void GenPrepDummy(Shr lambda) {
      GenIndShrsDummy(lambda);
      GenMultBatchDummy();
      GenIOBatchDummy();
    }

    void GenPrepDummy() {
      GenPrepDummy(Shr(0));
    }

    void SetNetwork(std::shared_ptr<scl::Network> network, std::size_t id) {
      mNetwork = network;
      mID = id;
      mParties = network->Size();
    }

    void SetThreshold(std::size_t threshold) {
      if ( mParties != threshold + 2*(mBatchSize - 1) + 1 )
	      throw std::invalid_argument("It must hold that n = t + 2(k-1) + 1");
      if ( mParties <= 2*threshold )
	      throw std::invalid_argument("It must hold that t < n/2");
      mThreshold = threshold;
    }

    // GENERATE F.I. PREPROCESSING

    void GenIndShrsPartiesSend();
    void GenIndShrsPartiesReceive();

    void GenUnpackedShrPartiesSend();
    void GenUnpackedShrPartiesReceive();

    void GenUnpackedMaskPartiesSend();
    void GenUnpackedMaskPartiesReceive();

    // Zero shares. Used for:
    // Inputs, Outputs, 3xMult
    void GenZeroPartiesSend();
    void GenZeroPartiesReceive();
    void GenZero2PartiesSend();
    void GenZero2PartiesReceive();

    // Zero shares for taking a product
    void GenZeroForProdPartiesSend();
    void GenZeroForProdPartiesReceive();

    // Execute the products
    void GenProdPartiesSendP1();
    void GenProdP1ReceivesAndSends();
    void GenProdPartiesReceive();

    // Mapping gates to preprocessed data
    void PopulateIndvShrs(std::shared_ptr<MultGate> gate) {
      if ( !gate->IsPadding() ) mMapIndShrs[gate] = mIndShrs[mCTRIndShrs++];
    }
    void PopulateIndvShrs(std::shared_ptr<InputGate> gate) {
      if ( !gate->IsPadding() ) mMapIndShrs[gate] = mIndShrs[mCTRIndShrs++];
    }
    void PopulateIndvShrs(std::shared_ptr<AddGate> gate) {
      // This is done in topological order, so previous (left&right) gates are already handled
      mMapIndShrs[gate] = mMapIndShrs[gate->GetLeft()] + mMapIndShrs[gate->GetRight()];
    }
    void PopulateIndvShrs(std::shared_ptr<OutputGate> gate) {
      mMapIndShrs[gate] = mMapIndShrs[gate->GetLeft()];
    }

    void PopulateMultBatches(std::shared_ptr<MultBatch> mult_batch) {
      mMapMultBatch[mult_batch] = mMultBatchFIPrep[mCTRMultBatches++];
    }

    void PopulateInputBatches(std::shared_ptr<InputBatch> input_batch) {
      mMapInputBatch[input_batch] = mIOBatchFIPrep[mCTRInOutBatches++];
    }

    void PopulateOutputBatches(std::shared_ptr<OutputBatch> output_batch) {
      mMapOutputBatch[output_batch] = mIOBatchFIPrep[mCTRInOutBatches++];
    }

    // Generate FD Prep from FI Prep
    
    // PREP INPUT BATCH
    void PrepInputPartiesSendOwner(std::shared_ptr<InputBatch> input_batch);
    
    void PrepInputOwnerReceives(std::shared_ptr<InputBatch> input_batch);


    // PREP OUTPUT BATCH
    void PrepOutputPartiesSendOwner(std::shared_ptr<OutputBatch> output_batch);
    
    void PrepOutputOwnerReceives(std::shared_ptr<OutputBatch> output_batch);

    // PREP MULT BATCH
    void PrepMultPartiesSendP1(std::shared_ptr<MultBatch> mult_batch);
    
    void PrepMultP1ReceivesAndSends();

    void PrepMultPartiesReceive(std::shared_ptr<MultBatch> mult_batch);


    // Populate shares of e_ij
    void PrecomputeEij() {
      vec_ZZ_pE e;
      e.SetLength(mBatch_l);
      vector<long> zero(mBatch_l, 0);
      rmfe.set_input(zero);
      ZZ_pE zeros = long2ZZpE(rmfe.get_result());

      for(std::size_t i = 0; i < mBatch_l; i++){
        vector<long> temp(mBatch_l, 0);
        temp[i] = 1;
        rmfe.set_input(temp);
        e[i] = long2ZZpE(rmfe.get_result());
      }

      for(std::size_t k = 0; k < mBatch_l; k++) {
      for (std::size_t i = 0; i < mBatch_m; i++) {
	      Shr shr = e[k];
	      for (std::size_t j = 0; j < mBatch_m; ++j) {
	        if (j == i) continue;
	        shr *= (Scheme_nm.alpha_set[mID+1] + Scheme_nm.alpha_set[j]) / (Scheme_nm.alpha_set[j] - Scheme_nm.alpha_set[i]); //TODO:
	      }
	      mSharesOfEij.emplace_back(shr);
      }
      }
    }

    // Populate shares of e_i
    void PrecomputeEi() {
      for (std::size_t i = 0; i < mBatch_m; i++) {
	    Shr shr(1);
	      for (std::size_t j = 0; j < mBatch_m; ++j) {
	        if (j == i) continue;
	        shr *= (Scheme_nm.alpha_set[mID+1] + Scheme_nm.alpha_set[j]) / (Scheme_nm.alpha_set[j] - Scheme_nm.alpha_set[i]);
	      }
	    mSharesOfEi.emplace_back(shr);
      }
    }

    // Populate vandermonde matrix
    void PrecomputeVandermonde() {
      mVandermonde.reserve(mParties);
      for (std::size_t i = 0; i < mParties; i++) {
	      mVandermonde.emplace_back(std::vector<Shr>());
	      mVandermonde[i].reserve(mThreshold + 1);
	      Shr entry(1);
	      for (std::size_t j = 0; j < mThreshold + 1; ++j) {
	        mVandermonde[i].emplace_back(entry);
	        entry *= gring.set_T[i];
	      }
      }
      mVandermonde_b = ZZpEmatrix2ZZpmatrix(mVandermonde, gr_degree, positive_integer);
    }

    // Maps
    std::map<std::shared_ptr<Gate>, Shr> mMapIndShrs;
    std::map<std::shared_ptr<MultBatch>, MultBatchFIPrep> mMapMultBatch;
    std::map<std::shared_ptr<InputBatch>, IOBatchFIPrep> mMapInputBatch;
    std::map<std::shared_ptr<OutputBatch>, IOBatchFIPrep> mMapOutputBatch;


    size_t GetIndShrsSize() const { return mIndShrs.size(); }
    
    size_t GetUnpackedShrsASize() const { return mUnpackedShrsA.size(); }
    size_t GetUnpackedShrsBSize() const { return mUnpackedShrsB.size(); }
    
    size_t GetUnpackedShrsABatchSize(size_t i) const { 
        return i < mUnpackedShrsA.size() ? mUnpackedShrsA[i].size() : 0; 
    }
    
    size_t GetUnpackedShrsBBatchSize(size_t i) const { 
        return i < mUnpackedShrsB.size() ? mUnpackedShrsB[i].size() : 0; 
    }
    
    size_t GetMultBatchFIPrepSize() const { return mMultBatchFIPrep.size(); }
    size_t GetIOBatchFIPrepSize() const { return mIOBatchFIPrep.size(); }

  private:
    // Sizes
    std::size_t mNIndShrs;
    std::size_t mNMultBatches;
    std::size_t mNInOutBatches;

    //packed scheme
    packed_shamir::scheme Scheme_nm;
    packed_shamir::scheme Scheme_n1;
    packed_shamir::scheme Scheme_m1;
    packed_shamir::scheme Scheme_t;

    // Counters
    std::size_t mCTRIndShrs;
    std::size_t mCTRMultBatches;
    std::size_t mCTRInOutBatches;

    std::size_t mBatchSize;
    std::size_t mBatch_l;
    std::size_t mBatch_m;


    // Shares of e_i for current party
    std::vector<Shr> mSharesOfEij; // len = batchsize, eij = i+j*l
    std::vector<Shr> mSharesOfEi; //len = batch_m
    std::vector<std::vector<FF>> mVandermonde_b;  // mParties*r x (mThreshold + 1)*r
    std::vector<std::vector<Shr>> mVandermonde; // mParties x (mThreshold + 1)
 

    // SHARINGS

    // Individual sharings
    std::vector<Shr> mIndShrs;

    // Multiplication batches
    std::vector<MultBatchFIPrep> mMultBatchFIPrep;

    // Input & Output batches
    std::vector<IOBatchFIPrep> mIOBatchFIPrep;

    // HELPERS TO GET THE F.I. PREP
    std::vector<std::vector<Shr>> mUnpackedShrsA; // idx = packed index
    std::vector<std::vector<Shr>> mUnpackedShrsB; // idx = packed index
    std::vector<std::vector<Shr>> mUnpackedShrsMask; // idx = packed index
    std::vector<std::vector<Shr>> mZeroProdShrs; // idx = packed index
    std::vector<std::vector<Shr>> mUnpackedShrsMaskPhiPsi; // idx = packed index

    

    // NETWORK-RELATED
    std::shared_ptr<scl::Network> mNetwork;
    std::size_t mID;
    std::size_t mParties;

    std::size_t mThreshold;

    scl::PRG mPRG;
  };

} // namespace tp

#endif  // CORRELATOR_H
