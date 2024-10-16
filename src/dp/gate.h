#ifndef GATE_H
#define GATE_H

#include <vector>
#include <assert.h>

#include "dp.h"

namespace dp {

  class Gate {
  public:

    // Individual share of the mask (lambda...lambda) associated to
    // the output of this gate.
    Shr GetShrLambda() { return mIndvShrLambdaC; };

    // Intended to be called by P1
    virtual FF GetMu() = 0;

    // Intended to be called by P1
    bool IsLearned() { return mLearned; };

    // For cleartext evaluation
    bool IsEvaluated() { return mEvaluated; }

    // For testing purposes. Sets mu to the given value
    void _SetDummyMu(FF mu) {
      mMu = mu;
      mLearned = true;
    }

    // To get indv_shr lambda 
    virtual Shr GetIndvShrLambda() = 0; 

    // Setting DN07 share
    //void SetDn07Share(FF share) {
    //  mDn07Share = share;
    //  mDn07Set = true;
    //}

    // To get share when run with DN07
    //virtual FF GetDn07Share() = 0; 

    // To get lambda when having fake preprocessing
    virtual FF GetDummyLambda() = 0; 
    
    // Get the cleartext value associated to the output of this gate,
    // if computed already (cleartext evaluation)
    virtual FF GetClear() = 0;

    std::shared_ptr<Gate> GetLeft() { return mLeft; }
    std::shared_ptr<Gate> GetRight() { return mRight; }

    bool IsPadding() { return mIsPadding; }

  protected:
    // Bool that indicates whether P1 learned Mu already
    bool mLearned = false;

    // Sharing of (lambda ... lambda)
    Shr mIndvShrLambdaC = Shr(0);
    bool mIndvShrLambdaCSet = false;

    // DN07-related
    //FF mDn07Share;
    //bool mDn07Set = false;

    // mu = value - lambda
    // Learned by P1 in the online phase
    FF mMu;

    // Actual lambda. Used for debugging purposes
    FF mLambda;
    bool mLambdaSet = false;

    // Parents
    std::shared_ptr<Gate> mLeft;
    std::shared_ptr<Gate> mRight;

    // For cleartext evaluation
    // Bool that indicates whether the gate has been evaluated
    bool mEvaluated = false;
    // Evaluation of the output wire of the gate
    FF mClear;

    bool mIsPadding = false;


    friend class MultBatch;
  };  

  class AddGate : public Gate {
  public:
    AddGate(std::shared_ptr<Gate> left, std::shared_ptr<Gate> right) {
      mLeft = left;
      mRight = right;
      mIndvShrLambdaC = left->GetShrLambda() + right->GetShrLambda();
    }

    FF GetMu() {
      if ( !mLearned ) {
	      mMu = mLeft->GetMu() + mRight->GetMu();
	      mLearned = true;
      }
      return mMu;
    }

    FF GetClear() {
      if ( !mEvaluated ) {
	      mClear = mLeft->GetClear() + mRight->GetClear();
	      mEvaluated = true;
      }
      return mClear;
    }

    Shr GetIndvShrLambda() {
      if ( !mIndvShrLambdaCSet ) {
	      mIndvShrLambdaC = mLeft->GetIndvShrLambda() + mRight->GetIndvShrLambda();
	      mIndvShrLambdaCSet = true;
      }
      return mIndvShrLambdaC;
    }    

//    FF GetDn07Share() {
//      if ( !mDn07Set ) {
//	mDn07Share = mLeft->GetDn07Share() + mRight->GetDn07Share();
//	mDn07Set = true;
//      }
//     return mDn07Share;
//    }    

    FF GetDummyLambda() {
      if ( !mLambdaSet ) {
	      mLambda = mLeft->GetDummyLambda() + mRight->GetDummyLambda();
	      mLambdaSet = true;
      }
      return mLambda;
    }    
  };
} // namespace tp

#endif  // GATE_H
