/**
 * \class StFemtoV0
 * \brief Keeps V0 daugheter IDs
 *
 * The class keeps indices of the positive and negative
 * daughters from StFemtoTrack StFemtoTrack array
 *
 * \author Egor Alpatov, Grigory Nigmatkulov; e-mail: nigmatkulov@gmail.com
 * \date November 20, 2018
 */

#ifndef StFemtoV0_h
#define StFemtoV0_h

// ROOT headers
#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "StRoot/KFParticle/KFParticle.h"

//________________
class StFemtoV0 : public TObject {

 public:
  /// Default costructor
  StFemtoV0();
  /// Copy constructor
  StFemtoV0(const StFemtoV0 &v0);
  /// Constructor that takes indices of positive and negative daughters
  /// in the StFemtoTrack array
  StFemtoV0(const Int_t& posId, const Int_t& negId, const bool& isParticle, const Float_t& decLength, const Float_t& dcaDaughter,
            const Float_t& primDCA, const Float_t& posPL, const Float_t& negPL, const TVector3& momentum, const Float_t& Minv)
  { mPosId = posId; mNegId = negId; mDecayLength = decLength; mDaughterDCA = dcaDaughter; mPrimaryDCA = primDCA; mPosPathLength = posPL;
    mNegPathLength = negPL; mMom = momentum; mInvMass = Minv; mIsParticle = isParticle;}
  /// Destructor
  virtual ~StFemtoV0()            { /* empty */ }
  /// Print V0 information
  void print();

  //
  // Setters
  //

  /// Set index of positive track in StFemtoTrack array
  void setPosId(const Int_t& pos) { mPosId = pos; }
  /// Set index of negative track in StFemtoTrack array
  void setNegId(const Int_t& neg) { mNegId = neg; }
  void setDecayLength(const Float_t& decLength) {mDecayLength = decLength; }
  void setDaughterDCA(const Float_t& dcaDaughter) {mDaughterDCA = dcaDaughter; }
  void setPrimaryDCA(const Float_t& primDCA) {mPrimaryDCA = primDCA; }
  void setPosPathLength(const Float_t& posPL) {mPosPathLength = posPL; }
  void setNegPathLength(const Float_t& negPL) {mNegPathLength = negPL; }
  void setMom(const TVector3& momentum) {mMom = momentum; }
  void setInvMass(const Float_t& Minv) {mInvMass = Minv; }
  void setIsParticle(const bool& isParticle) {mIsParticle = isParticle; }

  //
  // Getters
  //

  /// Return index of positive track in StFemtoTrack array
  Int_t posDaughterId() const     { return mPosId; }
  /// Return index of negative track in StFemtoTrack array
  Int_t negDaughterId() const     { return mNegId; }
  /// Return Decay Length
  Float_t decayLength() const     {return mDecayLength; }
  /// Return DCA between daughters
  Float_t daughterDCA() const	  {return mDaughterDCA; }
  /// Return DCA of V0 to PV
  Float_t primaryDCA() const	  {return mPrimaryDCA; }
  /// Return momentum of V0in decay point
  TVector3 Mom() const 		  {return mMom; }
  /// Return 4-Momentum of V0
  TLorentzVector FourMom(float Mass);
  /// Return pathLength for positive particle
  Float_t posPathLength() const	  {return mPosPathLength; }
  /// Return pathLength for negative particle
  Float_t negPathLength() const	  {return mNegPathLength; }
  ///Return invariant mass of daughter pair
  Float_t invMass() const 	  {return mInvMass; }
  //Return true if it is particle and false if antiparticle
  bool isParticle() const	  {return mIsParticle; }
  // Return PseudoKFParticle clone
  KFParticle kfparticle();   

 private:

  /// Index of positive track in StFemtoTrack array
  UShort_t mPosId;
  /// Index of negative track in StFemtoTrack array
  UShort_t mNegId;
  Float_t mDecayLength;
  Float_t mDaughterDCA;
  Float_t mPrimaryDCA;
  Float_t mPosPathLength;
  Float_t mNegPathLength;
  Float_t mInvMass;
  TVector3 mMom;
  bool mIsParticle;

  ClassDef(StFemtoV0, 1);
};

#endif // #define StFemtoV0_h
