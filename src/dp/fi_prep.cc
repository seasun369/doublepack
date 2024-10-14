#include "dp/correlator.h"

namespace dp {
  // GEN F.I. PREP
  void Correlator::GenIndShrsPartiesSend() {
    std::size_t degree = mParties - mBatch_m;
    std::size_t n_blocks = (mNIndShrs + (mThreshold + 1) -1) / (mThreshold + 1);
    for ( std::size_t block = 0; block < n_blocks; block++ ) {
      // 1 sample secret and shares
      //FF secret = random_ZZ_p();

      //vector<long> v(mBatch_l, conv<long>(secret));
      //rmfe.set_input(v);
      //ZZ_pE zeros = long2ZZpE(rmfe.get_result());

      for(std::size_t r = 0; r < gr_degree; r++){

      FF secret = random_ZZ_p();
      vec_ZZ_pE secrets;   //[s*1]_n-m
      secrets.SetLength(mBatch_m);  // 设置向量长度

      for (std::size_t i = 0; i < mBatch_m; i++) {
        secrets[i] = conv<Shr>(secret);  
      }
    
      //Vec secrets(std::vector<FF>(mBatchSize, secret));

      //auto poly = scl::details::EvPolyFromSecretsAndDegree(secrets, degree, mPRG);
      //auto shares = scl::details::SharesFromEvPoly(poly, mParties);

      vec_ZZ_pE shares = Scheme_nm.create_shares(secrets);
      //vector<vector<long>> shares_long;
      //for(int i=0; i<mParties; i++){
      //  ZzpE2Veclong( shares[i], shares_long[i], gr_degree);
      //}
      

      // 2 send shares
      for ( std::size_t party = 0; party < mParties; party++ ){
        mNetwork->Party(party)->Send(shares[party]);
      }
      }
    }
  }

  void Correlator::GenIndShrsPartiesReceive() {
    assert(mIndShrs.size() == 0);
    std::size_t n_blocks = (mNIndShrs + (mThreshold + 1) -1) / (mThreshold + 1);
    for ( std::size_t block = 0; block < n_blocks; block++ ) {
      // 1 receive shares
      std::vector<Shr> recv_shares;
      recv_shares.reserve(mParties*gr_degree);
      for(std::size_t r = 0; r < gr_degree; r++){
      for (std::size_t parties = 0; parties < mParties; parties++) {
	      Shr buffer;
        mNetwork->Party(parties)->Recv(buffer);
	      recv_shares.emplace_back(buffer);
      }
      }
      // 2 multiply by Vandermonde
      for ( std::size_t shr_idx = 0; shr_idx < (mThreshold+1)*gr_degree; shr_idx++ ){
	      Shr shr(0);
	      for ( std::size_t j = 0; j < mParties*gr_degree; j++ ){
	        shr += mVandermonde_b[j][shr_idx] * recv_shares[j];
	      }
	    mIndShrs.emplace_back(shr);
      }
    }
  }
/*
  generate Beaver triple
*/
  void Correlator::GenUnpackedShrPartiesSend() {
    std::size_t degree = mThreshold;
    std::size_t n_amount = 3*mNMultBatches; // 2 for the two factors, 1 for the multiplication
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);

    for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
      for ( std::size_t block = 0; block < n_blocks; block++ ) {
	      // 1 sample secret and shares
	      vector<long> aa(mBatch_l, 0);
        for(std::size_t i =0; i<mBatch_l; i++){
          FF a = random_ZZ_p();
          aa.emplace_back(conv<long>(a));
        }
        rmfe.set_input(aa);
        rmfe.RMFE_GR_PHI();
	      Shr secret = long2ZZpE(rmfe.get_result());
        //rmfe.RMFE_GR_PSI();

        vec_ZZ_pE shares = Scheme_t.create_one_shares(secret, pack_idx);

        //auto poly = scl::details::EvPolyFromSecretAndPointAndDegree(secret, FF(-pack_idx), degree, mPRG);
	      //auto shares = scl::details::SharesFromEvPoly(poly, mParties);
	
	      // 2 send shares
	      for ( std::size_t party = 0; party < mParties; party++ ){
	        mNetwork->Party(party)->Send(shares[party]);
	      }
      }
    }

  }
  void Correlator::GenUnpackedShrPartiesReceive() {
    std::size_t n_amount = 3*mNMultBatches;
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);
    mUnpackedShrsA.reserve(mBatch_m);
    mUnpackedShrsB.reserve(mBatch_m);
    //mUnpackedShrsMask.reserve(mBatch_m);
    //mUnpackedShrsMaskPhiPsi.reserve(mBatch_m);

    for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
      std::size_t ctr(0);
      mUnpackedShrsA.emplace_back(std::vector<Shr>());
      mUnpackedShrsA[pack_idx].reserve(mNMultBatches);
      mUnpackedShrsB.emplace_back(std::vector<Shr>());
      mUnpackedShrsB[pack_idx].reserve(mNMultBatches);
      //mUnpackedShrsMask.emplace_back(std::vector<Shr>());
      //mUnpackedShrsMask[pack_idx].reserve(mNMultBatches);
      //mUnpackedShrsMaskPhiPsi.emplace_back(std::vector<Shr>());
      //mUnpackedShrsMaskPhiPsi[pack_idx].reserve(mNMultBatches);

      for ( std::size_t block = 0; block < n_blocks; block++ ) {
	      // 1 receive shares
	      Vec recv_shares;
	      recv_shares.reserve(mParties);
	      for (std::size_t parties = 0; parties < mParties; parties++) {
	        Shr buffer;
	        mNetwork->Party(parties)->Recv(buffer);
	        recv_shares.emplace_back(buffer);
	      }
	      // 2 multiply by Vandermonde
	      for ( std::size_t shr_idx = 0; shr_idx < mThreshold+1; shr_idx++ ){
	        Shr shr(0);
	        for ( std::size_t j = 0; j < mParties; j++ ){
	          shr += mVandermonde[j][shr_idx] * recv_shares[j];
	        }
	        if (ctr < mNMultBatches) mUnpackedShrsA[pack_idx].emplace_back(shr);
	        if ( (mNMultBatches <= ctr) && (ctr < 2*mNMultBatches) ) mUnpackedShrsB[pack_idx].emplace_back(shr);
	        //if ( (2*mNMultBatches <= ctr) && (ctr < 3*mNMultBatches) ) mUnpackedShrsMask[pack_idx].emplace_back(shr);
	        ctr++;
	      }
      }
    }
    // Create Mult-related data
    for ( std::size_t i = 0; i < mNMultBatches; i++ ) {
      MultBatchFIPrep data(Shr(0));
      for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
	      data.mShrA += mSharesOfEi[pack_idx] * mUnpackedShrsA[pack_idx][i];
	      data.mShrB += mSharesOfEi[pack_idx] * mUnpackedShrsB[pack_idx][i];
      }
      mMultBatchFIPrep.emplace_back(data);
    }
  }

/*
  generate r and phi(psi(r))
*/
  void Correlator::GenUnpackedMaskPartiesSend() {
    std::size_t degree = mThreshold;
    std::size_t n_amount = 3*mNMultBatches; // 2 for the two factors, 1 for the multiplication
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);

    for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
      for ( std::size_t block = 0; block < n_blocks; block++ ) {
	      // 1 sample secret and shares
	      //vector<long> aa(mBatch_l, 0);
        //for(std::size_t i =0; i<mBatch_l; i++){
        //  FF a = random_ZZ_p();
        //  aa.emplace_back(conv<long>(a));
        //}
        //Shr a = random_ZZ_pE();

        //rmfe.set_input(aa);
        //rmfe.RMFE_GR_PHI();
        //vector<long> result1 = rmfe.get_result();
	      Shr secret1 =  random_ZZ_pE();
        vector<long> result1;
        ZzpE2Veclong(secret1, result1, s);
        rmfe.RMFE_GR_PSI(result1);
        rmfe.RMFE_GR_PHI();
        vector<long> result2 = rmfe.get_result();
        Shr secret2 = long2ZZpE(result2);

        vec_ZZ_pE shares1 = Scheme_t.create_one_shares(secret1, pack_idx);
        vec_ZZ_pE shares2 = Scheme_t.create_one_shares(secret2, pack_idx);

        //auto poly = scl::details::EvPolyFromSecretAndPointAndDegree(secret, FF(-pack_idx), degree, mPRG);
	      //auto shares = scl::details::SharesFromEvPoly(poly, mParties);
	
	      // 2 send shares
	      for ( std::size_t party = 0; party < mParties; party++ ){
	        mNetwork->Party(party)->Send(shares1[party]);
          mNetwork->Party(party)->Send(shares2[party]);
	      }
      }
    }

  }
  void Correlator::GenUnpackedMaskPartiesReceive() {
    std::size_t n_amount = 3*mNMultBatches;
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);
    //mUnpackedShrsA.reserve(mBatch_m);
    //mUnpackedShrsB.reserve(mBatch_m);
    mUnpackedShrsMask.reserve(mBatch_m);
    mUnpackedShrsMaskPhiPsi.reserve(mBatch_m);

    for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
      std::size_t ctr(0);
      //mUnpackedShrsA.emplace_back(std::vector<Shr>());
      //mUnpackedShrsA[pack_idx].reserve(mNMultBatches);
      //mUnpackedShrsB.emplace_back(std::vector<Shr>());
      //mUnpackedShrsB[pack_idx].reserve(mNMultBatches);
      mUnpackedShrsMask.emplace_back(std::vector<Shr>());
      mUnpackedShrsMask[pack_idx].reserve(mNMultBatches);
      mUnpackedShrsMaskPhiPsi.emplace_back(std::vector<Shr>());
      mUnpackedShrsMaskPhiPsi[pack_idx].reserve(mNMultBatches);

      for ( std::size_t block = 0; block < n_blocks; block++ ) {
	      // 1 receive shares
	      Vec recv_shares1,recv_shares2;
	      recv_shares1.reserve(mParties);
        recv_shares2.reserve(mParties);
	      for (std::size_t parties = 0; parties < mParties; parties++) {
	        Shr buffer1,buffer2;
	        mNetwork->Party(parties)->Recv(buffer1);
          mNetwork->Party(parties)->Recv(buffer2);
	        recv_shares1.emplace_back(buffer1);
          recv_shares2.emplace_back(buffer1);
	      }
	      // 2 multiply by Vandermonde
	      for ( std::size_t shr_idx = 0; shr_idx < mThreshold+1; shr_idx++ ){
	        Shr shr1(0),shr2(0);
	        for ( std::size_t j = 0; j < mParties; j++ ){
	          shr1 += mVandermonde[j][shr_idx] * recv_shares1[j];
	        }
          for ( std::size_t j = 0; j < mParties; j++ ){
	          shr2 += mVandermonde[j][shr_idx] * recv_shares2[j];
	        }
	        if (ctr < mNMultBatches) {mUnpackedShrsMask[pack_idx].emplace_back(shr1);mUnpackedShrsMaskPhiPsi[pack_idx].emplace_back(shr2);}
	        //if ( (mNMultBatches <= ctr) && (ctr < 2*mNMultBatches) ) mUnpackedShrsB[pack_idx].emplace_back(shr);
	        //if ( (2*mNMultBatches <= ctr) && (ctr < 3*mNMultBatches) ) mUnpackedShrsMask[pack_idx].emplace_back(shr);
	        ctr++;
	      }
      }
    }
    // Create Mult-related data
    //for ( std::size_t i = 0; i < mNMultBatches; i++ ) {
    //  MultBatchFIPrep data(Shr(0));
    //  for ( std::size_t pack_idx = 0; pack_idx < mBatchSize; pack_idx++ ) {
	  //    data.mShrA += mSharesOfEi[pack_idx] * mUnpackedShrsA[pack_idx][i];
	  //    data.mShrB += mSharesOfEi[pack_idx] * mUnpackedShrsB[pack_idx][i];
    //  }
    //  mMultBatchFIPrep.emplace_back(data);
    //}
  }

  // Zero shares. Used for:
  // Inputs, Outputs, 3xMult but only [0]
  void Correlator::GenZeroPartiesSend() {
    std::size_t degree = mParties - 1;
    std::size_t n_amount = 3*mNMultBatches + mNInOutBatches;
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);
    for ( std::size_t block = 0; block < n_blocks; block++ ) {
      // 1 sample secret and shares
      vec_ZZ_pE secrets;
      secrets.SetLength(mBatch_m);  // 设置向量长度

      // 将每个元素初始化为零
      for (long i = 0; i < mBatch_m; i++) {
        secrets[i] = ZZ_pE::zero();  // 初始化为全零元素
      }


      //auto poly = scl::details::EvPolyFromSecretsAndDegree(secrets, degree, mPRG);
      //auto shares = scl::details::SharesFromEvPoly(poly, mParties);
      vec_ZZ_pE shares = Scheme_n1.create_shares(secrets);
	
      // 2 send shares
      for ( std::size_t party = 0; party < mParties; party++ ){
	      mNetwork->Party(party)->Send(shares[party]);
      }
    }
  }

  void Correlator::GenZeroPartiesReceive() {
    std::size_t n_amount = 3*mNMultBatches + mNInOutBatches;
    std::size_t ctr(0);
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);
    for ( std::size_t block = 0; block < n_blocks; block++ ) {
      // 1 receive shares
      Vec recv_shares;
      recv_shares.reserve(mParties);
      for (std::size_t parties = 0; parties < mParties; parties++) {
	      Shr buffer;
	      mNetwork->Party(parties)->Recv(buffer);
	      recv_shares.emplace_back(buffer);
      }
      // 2 multiply by Vandermonde
      for ( std::size_t shr_idx = 0; shr_idx < mThreshold+1; shr_idx++ ){
	      Shr shr(0);
	      for ( std::size_t j = 0; j < mParties; j++ ){
	        shr += mVandermonde[j][shr_idx] * recv_shares[j];
	      }
	      if (ctr < mNMultBatches) mMultBatchFIPrep[ctr].mShrO1 = shr; 
	      if ((mNMultBatches <= ctr) && (ctr < 2*mNMultBatches)) mMultBatchFIPrep[ctr - mNMultBatches].mShrO2 = shr;
        //else continue;
	      //if ((2*mNMultBatches <= ctr) && (ctr < 3*mNMultBatches)) mMultBatchFIPrep[ctr - 2*mNMultBatches].mShrO3 = shr;
	      //else {
	      //  IOBatchFIPrep tmp;
	      //  tmp.mShrO = shr;
	      //  mIOBatchFIPrep.emplace_back(tmp);   // should be <0>n-1
	      //}
	    ctr++;
      }
    }
  }

  // Zero shares. Used for:
  // Inputs, Outputs, 3xMult but for <0>
  void Correlator::GenZero2PartiesSend() {
    std::size_t degree = mParties - 1;
    std::size_t n_amount = 3*mNMultBatches + mNInOutBatches;
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);
    for ( std::size_t block = 0; block < n_blocks; block++ ) {
      // 1 sample secret and shares
      vec_ZZ_pE kernel;
      kernel.SetLength(mBatch_m);

      rmfe.get_phi_kernel(kernel, mBatch_m);

      vec_ZZ_pE secrets;
      secrets.SetLength(mBatch_m);
      for(long i=0; i<mBatch_m; i++){
        secrets[i] = kernel[i]; //should sample from kernel
      }

      //auto poly = scl::details::EvPolyFromSecretsAndDegree(secrets, degree, mPRG);
      //auto shares = scl::details::SharesFromEvPoly(poly, mParties);
      vec_ZZ_pE shares = Scheme_n1.create_shares(secrets);
	
      // 2 send shares
      for ( std::size_t party = 0; party < mParties; party++ ){
	      mNetwork->Party(party)->Send(shares[party]);
      }
    }
  }

  void Correlator::GenZero2PartiesReceive() {
    std::size_t n_amount = 3*mNMultBatches + mNInOutBatches;
    std::size_t ctr(0);
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);
    for ( std::size_t block = 0; block < n_blocks; block++ ) {
      // 1 receive shares
      Vec recv_shares;
      recv_shares.reserve(mParties);
      for (std::size_t parties = 0; parties < mParties; parties++) {
	      Shr buffer;
	      mNetwork->Party(parties)->Recv(buffer);
	      recv_shares.emplace_back(buffer);
      }
      // 2 multiply by Vandermonde
      for ( std::size_t shr_idx = 0; shr_idx < mThreshold+1; shr_idx++ ){
	      Shr shr(0);
	      for ( std::size_t j = 0; j < mParties; j++ ){
	        shr += mVandermonde[j][shr_idx] * recv_shares[j];
	      }
	      if (ctr < mNMultBatches) mMultBatchFIPrep[ctr].mShrO3 = shr; 
	      if ((mNMultBatches <= ctr) && (ctr < mNMultBatches+mNInOutBatches)) {
          IOBatchFIPrep tmp;
	        tmp.mShrO = shr;
	        mIOBatchFIPrep.emplace_back(tmp);
        }
	    ctr++;
      }
    }
  }

  void Correlator::GenZeroForProdPartiesSend() {
    std::size_t degree = mParties - 1;
    std::size_t n_amount = mNMultBatches; 
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);

    for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
      for ( std::size_t block = 0; block < n_blocks; block++ ) {
	    // 1 sample secret and shares
	    vec_ZZ_pE kernel;
      kernel.SetLength(mBatch_m);

      rmfe.get_phi_kernel(kernel, mBatch_m);

      vec_ZZ_pE secrets;
      secrets.SetLength(mBatch_m);
      for(long i=0; i<mBatch_m; i++){
        secrets[i] = kernel[i]; //should sample from kernel
      }

      //auto poly = scl::details::EvPolyFromSecretsAndDegree(secrets, degree, mPRG);
      //auto shares = scl::details::SharesFromEvPoly(poly, mParties);
      vec_ZZ_pE shares = Scheme_n1.create_shares(secrets);
	
	      // 2 send shares
	      for ( std::size_t party = 0; party < mParties; party++ ){
	        mNetwork->Party(party)->Send(shares[party]);
	      }
      }
    }
  }

  void Correlator::GenZeroForProdPartiesReceive() {
    std::size_t n_amount = mNMultBatches;
    std::size_t n_blocks = (n_amount + (mThreshold + 1) -1) / (mThreshold + 1);
    mZeroProdShrs.reserve(mBatch_m);

    for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
      mZeroProdShrs.emplace_back(std::vector<Shr>());
      mZeroProdShrs[pack_idx].reserve(mNMultBatches);

      for ( std::size_t block = 0; block < n_blocks; block++ ) {
	      // 1 receive shares
	      Vec recv_shares;
	      recv_shares.reserve(mParties);
	      for (std::size_t parties = 0; parties < mParties; parties++) {
	        Shr buffer;
	        mNetwork->Party(parties)->Recv(buffer);
	        recv_shares.emplace_back(buffer);
	      }
	      // 2 multiply by Vandermonde
	      for ( std::size_t shr_idx = 0; shr_idx < mThreshold+1; shr_idx++ ){
	        Shr shr(0);
	        for ( std::size_t j = 0; j < mParties; j++ ){
	          shr += mVandermonde[j][shr_idx] * recv_shares[j];
	        }
	        mZeroProdShrs[pack_idx].emplace_back(shr);
	      }
      }
    }
  }

  void Correlator::GenProdPartiesSendP1() {
    for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
      for ( std::size_t batch = 0; batch < mNMultBatches; batch++ ) {
        // 1. Gather shares
	      Shr share = mUnpackedShrsA[pack_idx][batch] * mUnpackedShrsB[pack_idx][batch]\
	      + mUnpackedShrsMask[pack_idx][batch] + mZeroProdShrs[pack_idx][batch];

	      // 2. send shares
	      mNetwork->Party(0)->Send(share);
      }
    }
  }

  void Correlator::GenProdP1ReceivesAndSends() {
    if ( mID == 0 ) {
      for ( std::size_t pack_idx = 0; pack_idx < mBatch_m; pack_idx++ ) {
	      for ( std::size_t batch = 0; batch < mNMultBatches; batch++ ) {
	      // 1. Receive shares
	        Vec recv_shares;
	        recv_shares.reserve(mParties);
	        for (std::size_t parties = 0; parties < mParties; parties++) {
	          Shr buffer;
	          mNetwork->Party(parties)->Recv(buffer);
	          recv_shares.emplace_back(buffer);
	        }

	        // 2. Reconstruct
	        //auto secret = SecretFromPointAndShares(FF(-pack_idx), recv_shares);
          vec_ZZ_pE rec;
          rec.SetLength(mParties);
          for(std::size_t parties = 0; parties < mParties; parties++){
            rec[parties] = recv_shares[parties];
          }

          Shr secret = Scheme_n1.reconstruct_one_shares(rec, pack_idx);

	        // 3. Send back (w. optimization of zero-shares)
	        Vec y_points;
	        y_points.reserve(mThreshold+1);
	        y_points.emplace_back(secret);
	        for (std::size_t i = 1; i < mThreshold+1; ++i) y_points.emplace_back(Shr(0));

	        Vec x_points;
	        x_points.reserve(mThreshold+1);
	        x_points.emplace_back(Scheme_n1.beta_set(pack_idx));
	        for (std::size_t i = 1; i < mThreshold+1; ++i) x_points.emplace_back(Scheme_n1.beta_set(i+mBatch_m));

	        //auto poly = scl::details::EvPolynomial<FF>(x_points, y_points);
	        //auto shares_to_send = scl::details::SharesFromEvPoly(poly, mParties);
          vec_ZZ_pE shares = Scheme_n1.create_shares_with_points(x_points,y_points);


	        for ( std::size_t i = mThreshold; i < mParties; i++ ) {
	          mNetwork->Party(i)->Send(shares[i]);
	        }
	      }
      }
    }
  }

  void Correlator::GenProdPartiesReceive() {
    std::vector<std::vector<Shr>> shares_prod;
    shares_prod.reserve(mBatchSize);

    for ( std::size_t pack_idx = 0; pack_idx < mBatchSize; pack_idx++ ) {
      shares_prod.emplace_back(std::vector<Shr>());
      shares_prod[pack_idx].reserve(mNMultBatches);
      for ( std::size_t batch = 0; batch < mNMultBatches; batch++ ) {
	      Shr recv;
	      if (mID >= mThreshold) {
	      // 1. Receive secret
	        mNetwork->Party(0)->Recv(recv);
	      } else {
	        recv = Shr(0);
	      }
	       
	      // 2. Compute shares
	      Shr shr_prod = recv - mUnpackedShrsMaskPhiPsi[pack_idx][batch];
	      shares_prod[pack_idx].emplace_back(shr_prod);
      }
    }

    for ( std::size_t i = 0; i < mNMultBatches; i++ ) {
      for ( std::size_t pack_idx = 0; pack_idx < mBatchSize; pack_idx++ ) {
	      mMultBatchFIPrep[i].mShrC += mSharesOfEi[pack_idx] * shares_prod[pack_idx][i];
      }
    }
  }
}
