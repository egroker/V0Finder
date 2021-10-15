#ifndef StPicoV0Finder_h
#define StPicoV0Finder_h

#include <string>

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"

#include "StFemtoV0/StFemtoV0.h"
#include "StFemtoV0/StFemtoXi.h"

class StPicoV0Finder : public TObject {

 public:
 StFemtoV0 V0Builder(int protonId, int pionId, bool isParticle);
 StFemtoXi XiBuilder(int partId, StFemtoV0 v0, bool isParticle);
 StPicoV0Finder(StPicoDst *wholeDst, StPicoEvent *wholeEvent);
 virtual ~StPicoV0Finder() {}
 void SetDCAtoPVCut(float dca) {mDCAtoPVCut = dca;}
 void SetDaughterDCACut(float daughDca) {mDaughterDCACut = daughDca;}
 void SetDecayLengthCut(float decLength) {mDecayLengthCut = decLength;}

 void SetXiDCAtoPVTopCut(float dcaTop) {mXiDCAtoPVTopCut = dcaTop;}
 void SetXiDCAtoPVLowCut(float dcaLow) {mXiDCAtoPVLowCut = dcaLow;}
 void SetXiDaughterDCACut(float daughDca) {mXiDaughterDCACut = daughDca;}
 void SetXiDecayLengthCut(float decLength) {mXiDecayLengthCut = decLength;}
 void SetXiMassWindow(float massWindow) {mXiMassWindow = massWindow;}
 StPicoPhysicalHelix helix(StFemtoV0 v0, float field);
 TVector3 decayPoint(StFemtoV0 v0, TVector3 primvertex);
 float dcaToPV() {return mLamDca;}

 private:
 StPicoDst *dst;
 StPicoEvent *event;
 
 float mDCAtoPVCut, mDaughterDCACut, mDecayLengthCut;
 float mXiDCAtoPVTopCut, mXiDCAtoPVLowCut, mXiDaughterDCACut, mXiDecayLengthCut;
 float mXiMassWindow;

 bool PionCut(StPicoTrack *pion);
 bool ProtonCut(StPicoTrack *proton);
 bool LambdaCut(float DCAtoPV, float DaughterDCA, float DecayLength, float PointAway, float MomPerp, float LamMass);
 bool XiCut(float DCAtoPV, float DaughterDCA, float DecayLength, float PointAway, float MomPerp, float XiMass);
 bool XiDaughterCut(float DCAtoPV, float DaughterDCA, float DecayLength, float PointAway, float MomPerp, float XiMass);

 const float MassProton = 0.9382720;
 const float MassPion = 0.139570;
 const float MassKaon = 0.4936770;
 const float MassElectron = 0.5110e-03;
 const float MassD = 1.87561339;
 const float MassT = 2.80925;
 const float MassHe3 = 2.80923;
 const float MassLambda = 1.115683;
 const float MassXi = 1.32171;
 float mLamDca, mLamDaughterDCA, mLamDecLength, mLamMass, mLamPointAway, mLamMom;
 TVector3 mLamStart;
 ClassDef(StPicoV0Finder, 1);
};

#endif
