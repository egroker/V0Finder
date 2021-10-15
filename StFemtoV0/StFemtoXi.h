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

#ifndef StFemtoXi_h
#define StFemtoXi_h

// ROOT headers
#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "StFemtoV0.h"
//________________
class StFemtoXi : public TObject {

 public:
  /// Default costructor
  StFemtoXi();
  /// Copy constructor
  StFemtoXi(const StFemtoXi &xi);
  /// Constructor that takes indices of positive and negative daughters
  /// in the StFemtoTrack array
  StFemtoXi(const Int_t& partId, StFemtoV0 &v0, const bool& isParticle, const Float_t& decLength, const Float_t& dcaDaughter,
            const Float_t& primDCA, const Float_t& posPL, const Float_t& negPL, const TVector3& momentum, const Float_t& Minv)
  { mPartId = partId; mV0 = v0; mDecayLength = decLength; mDaughterDCA = dcaDaughter; mPrimaryDCA = primDCA; mPosPathLength = posPL;
    mNegPathLength = negPL; mMom = momentum; mInvMass = Minv; mIsParticle = isParticle;}
  /// Destructor
  virtual ~StFemtoXi()            { /* empty */ }
  /// Print V0 information
  void print();

  //
  // Setters
  //

  /// Set index of positive track in StFemtoTrack array
  void setPartId(const Int_t& id) { mPartId = id; }
  /// Set index of negative track in StFemtoTrack array
  void setV0(StFemtoV0 &v0) { mV0 = v0; }
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
  Int_t DaughterId() const     { return mPartId; }
  /// Return index of negative track in StFemtoTrack array
  StFemtoV0 getV0() const     { return mV0; }
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
  bool isParticle() const	  {return mIsParticle; }

 private:

  /// Index of positive track in StFemtoTrack array
  UShort_t mPartId;
  /// Index of negative track in StFemtoTrack array
  StFemtoV0 mV0;
  Float_t mDecayLength;
  Float_t mDaughterDCA;
  Float_t mPrimaryDCA;
  Float_t mPosPathLength;
  Float_t mNegPathLength;
  Float_t mInvMass;
  TVector3 mMom;
  bool mIsParticle;

  ClassDef(StFemtoXi, 1);
};

#endif // #define StFemtoXi_h
