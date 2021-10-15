//
// The StFemtoV0 class holds indices of positive and
// negative daughters in StFemtoTrack array
//

// C++ headers
#include <iostream>

// FemtoDst headers
#include "StFemtoXi.h"
#include "StFemtoV0.h"

// ROOT headers
#include "TLorentzVector.h"
ClassImp(StFemtoXi)

//________________
StFemtoXi::StFemtoXi() : TObject(), mPartId(0), mDecayLength(0), mDaughterDCA(0),
                         mPrimaryDCA(0), mPosPathLength(0), mNegPathLength(0), mInvMass(0), mIsParticle(0) {
  /* empty */
}

//________________
StFemtoXi::StFemtoXi(const StFemtoXi &v0) :
  TObject(), mPartId(v0.mPartId), mV0(v0.mV0), mDecayLength(v0.mDecayLength), mDaughterDCA(v0.mDaughterDCA),
                         mPrimaryDCA(v0.mPrimaryDCA), mPosPathLength(v0.mPosPathLength), mNegPathLength(v0.mNegPathLength),
                         mInvMass(v0.mInvMass), mMom(v0.mMom), mIsParticle(v0.mIsParticle) {
  /* empty */
}

//________________
void StFemtoXi::print() {
  std::cout << "V0 daughters. PosId: " << mPartId << std::endl;
  std::cout << "Minv: " << mInvMass << "DCA to PV: " << mPrimaryDCA << std::endl;
  std::cout << "PL1 : " << mPosPathLength << "PL2 : " << mNegPathLength << std::endl;  
}

//______________
TLorentzVector StFemtoXi::FourMom(float Mass) {
  TLorentzVector FourMom(mMom.X(), mMom.Y(), mMom.Z(), sqrt(mMom.Mag2() + Mass*Mass));
  return FourMom;
}

