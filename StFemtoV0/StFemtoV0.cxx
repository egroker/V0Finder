//
// The StFemtoV0 class holds indices of positive and
// negative daughters in StFemtoTrack array
//

// C++ headers
#include <iostream>

// FemtoDst headers
#include "StFemtoV0.h"

// ROOT headers
#include "TLorentzVector.h"
ClassImp(StFemtoV0)

//________________
StFemtoV0::StFemtoV0() : TObject(), mPosId(0), mNegId(0), mDecayLength(0), mDaughterDCA(0),
                         mPrimaryDCA(0), mPosPathLength(0), mNegPathLength(0), mInvMass(0), mIsParticle(0) {
  /* empty */
}

//________________
StFemtoV0::StFemtoV0(const StFemtoV0 &v0) :
  TObject(), mPosId(v0.mPosId), mNegId(v0.mNegId), mDecayLength(v0.mDecayLength), mDaughterDCA(v0.mDaughterDCA),
                         mPrimaryDCA(v0.mPrimaryDCA), mPosPathLength(v0.mPosPathLength), mNegPathLength(v0.mNegPathLength),
                         mInvMass(v0.mInvMass), mMom(v0.mMom), mIsParticle(v0.mIsParticle) {
  /* empty */
}

//________________
void StFemtoV0::print() {
  std::cout << "V0 daughters. PosId: " << mPosId << " NegId: "
            << mNegId << std::endl;
  std::cout << "Minv: " << mInvMass << "DCA to PV: " << mPrimaryDCA << std::endl;
  std::cout << "PL1 : " << mPosPathLength << "PL2 : " << mNegPathLength << std::endl;  
}

//______________
TLorentzVector StFemtoV0::FourMom(float Mass) {
  TLorentzVector FourMom(mMom.X(), mMom.Y(), mMom.Z(), sqrt(mMom.Mag()*mMom.Mag() + Mass*Mass));
  return FourMom;
}

//________________
