#include "dp/mult_gate.h"

namespace dp {
  void MultBatch::P1Sends() {
    if ( mID == 0 ) {

      // 1. P1 assembles mu_A and mu_B

      Vec mu_alpha;
      Vec mu_beta;
      // Here is where the permutation happens!
      for (auto gate : mMultGatesPtrs) {
	      mu_alpha.emplace_back(gate->GetLeft()->GetMu());
	      mu_beta.emplace_back(gate->GetRight()->GetMu());
      }
      
      // 2. P1 generates shares of mu_A and mu_B
      
      //auto poly_A = scl::details::EvPolyFromSecretsAndDegree(mu_alpha, mBatchSize-1, mPRG);
      //Vec shares_A = scl::details::SharesFromEvPoly(poly_A, mParties);
      vec_ZZ_pE shares_A = Scheme_m1.create_shares(mu_alpha);

      //auto poly_B = scl::details::EvPolyFromSecretsAndDegree(mu_beta, mBatchSize-1, mPRG);
      //Vec shares_B = scl::details::SharesFromEvPoly(poly_B, mParties);
      vec_ZZ_pE shares_B = Scheme_m1.create_shares(mu_beta);

      // 3. P1 sends the shares

      for (std::size_t i = 0; i < mParties; ++i) {
	      mNetwork->Party(i)->Send(shares_A[i]);
	      mNetwork->Party(i)->Send(shares_B[i]);
      }
    }
  }

  void MultBatch::PartiesReceive() {
      Shr shr_mu_A;
      Shr shr_mu_B;

      mNetwork->Party(0)->Recv(shr_mu_A);
      mNetwork->Party(0)->Recv(shr_mu_B);

      mPackedShrMuA = shr_mu_A;
      mPackedShrMuB = shr_mu_B;
    }

  void MultBatch::PartiesSend() {
      // Compute share
      Shr shr_mu_C;
      shr_mu_C = mPackedShrMuB * mPackedShrLambdaA + mPackedShrMuA * mPackedShrLambdaB + \
	    mPackedShrMuA * mPackedShrMuB + mPackedShrDeltaC;

      // Send to P1
      mNetwork->Party(0)->Send(shr_mu_C); 
    }

  void MultBatch::P1Receives() {
    if (mID == 0) {
	    vec_ZZ_pE mu_gamma;
      Vec shares;
	      for (std::size_t i = 0; i < mParties; ++i) {
        FF shr_mu_C;
        mNetwork->Party(i)->Recv(shr_mu_C);
        shares.emplace_back(shr_mu_C);
      }
	    //mu_gamma = scl::details::SecretsFromSharesAndLength(shares, mBatchSize);
      mu_gamma = Scheme_m1.packed_reconstruct_shares(shares); //TODO: use rmfe: psi^prime

	  // P1 updates the mu for the gates in the current batch
	      for (std::size_t i = 0; i < mBatchSize; i++) {
	        mMultGatesPtrs[i]->mMu = mu_gamma[i];
	        mMultGatesPtrs[i]->mLearned = true;
	      }
      }
    }
    
} // namespace tp
  
