#include <vector>

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"

#include "StFemtoV0/StFemtoV0.h"
#include "StFemtoV0/StFemtoXi.h"

#include "StPicoV0Finder.h"
#include "TLorentzVector.h"
StPicoV0Finder::StPicoV0Finder(StPicoDst *wholeDst, StPicoEvent *wholeEvent):
  dst( wholeDst ), event( wholeEvent ) {
  mDCAtoPVCut = 10.; mDaughterDCACut = 2.; mDecayLengthCut = 2.;
  mXiDCAtoPVTopCut = 20.; mXiDCAtoPVLowCut = 2.; mDaughterDCACut = 2.; mDecayLengthCut = 2.;
  /* empty */
  }

StFemtoV0 StPicoV0Finder::V0Builder(int protonId, int pionId, bool isParticle) {
  StFemtoV0 vzero;
  float DCAtoPV1, DCAtoPV2;  //dca vector to primary vertex
  float px1, px2, px3, py1, py2, py3, pz1, pz2, pz3, P1, P2, P3, E1, E2, E3;
  std::pair<double,double> lengths;
  TVector3 closest1, closest2, LambdaDecMom;
  double particleDCA;
  TVector2 KaonPhiVec;
  float KaonPhi;
  TVector3 v0PointDecay, MomClosest1, MomClosest2, PointLambdaToVtx, mDca;
  double LambdaDCA;
  float Minv, LambdaMass, KinkMass, AntiKinkMass;
  mLamDca = 0;        
        StPicoTrack *mGlobTrack;
        mGlobTrack = dst->track(protonId); //Getting Positive Particle
        if(!mGlobTrack) return vzero;
        if(!ProtonCut(mGlobTrack)) return vzero;
        px1 = mGlobTrack->gMom().x();
        py1 = mGlobTrack->gMom().y();
        pz1 = mGlobTrack->gMom().z();
        TVector3 pvec1(px1,py1,pz1);
        TVector3 origin1(mGlobTrack->origin().x(),
                                      mGlobTrack->origin().y(),
                                      mGlobTrack->origin().z());
        P1 = mGlobTrack->gMom().Mag();
        E1 = sqrt(P1*P1 + MassProton*MassProton);
        StPicoPhysicalHelix helix1(pvec1, origin1, event->bField()/pow(10,14) , mGlobTrack->charge());
        DCAtoPV1 = mGlobTrack->gDCA(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());

        TLorentzVector PartPos1(px1, py1, pz1, E1);

          mGlobTrack = dst->track(pionId ); //Getting Negative Particle
          if(!mGlobTrack) return vzero;
	  if (!PionCut(mGlobTrack)) return vzero;
          px2 = mGlobTrack->gMom().x();
          py2 = mGlobTrack->gMom().y();
          pz2 = mGlobTrack->gMom().z();
          P2 = mGlobTrack->gMom().Mag();
          E2 = sqrt(P2*P2 + MassPion*MassPion);
          TVector3 pvec2(px2,py2,pz2);
          TVector3 origin2(mGlobTrack->origin().x(),
                           mGlobTrack->origin().y(),
                           mGlobTrack->origin().z());
          StPicoPhysicalHelix helix2(pvec2,origin2,event->bField()/pow(10,14), mGlobTrack->charge());
          DCAtoPV2 = mGlobTrack->gDCA(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
          TLorentzVector PartNeg1(px2, py2, pz2, E2);

          TVector3 closestDiff12;
          lengths = helix1.pathLengths(helix2);         //dca of one helix to the other
          closest1 = helix1.at(lengths.first);
          closest2 = helix2.at(lengths.second);
          closestDiff12 = closest2 - closest1;
          particleDCA = closestDiff12.Mag();

          v0PointDecay = (closest1 + closest2)*0.5; // Calculating DCA to PV
          MomClosest1 = helix1.momentumAt(lengths.first, event->bField()/pow(10,14));
          MomClosest2 = helix2.momentumAt(lengths.second, event->bField()/pow(10,14));
          LambdaDecMom = (MomClosest1 + MomClosest2);

          TVector3 backdir = LambdaDecMom*(1/(LambdaDecMom.Mag()));
          TVector3 primvertex(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
          TVector3 LambdaPrimDiff = v0PointDecay - primvertex;

          double MinFac = - ((LambdaPrimDiff.x())*(backdir.x())
                         + (LambdaPrimDiff.y())*(backdir.y())
                         + (LambdaPrimDiff.z())*(backdir.z()))/(backdir.Mag()); 
          float angle = (backdir*MinFac).Angle(LambdaPrimDiff);
          float cosinus = TMath::Cos(angle);
          float sinus = TMath::Sin(angle);
          float DLcos = (LambdaPrimDiff*cosinus).Mag();
          PointLambdaToVtx = v0PointDecay + MinFac*backdir;
          TVector3 VectLambdaDCA = PointLambdaToVtx - primvertex;
          LambdaDCA = VectLambdaDCA.Mag(); //(k0PointDecay-primvertex).Mag()*sinus; //(KaonPrimDiff*sinus).Mag(); 
	  mLamDca = LambdaDCA; 
          float LambdaDCAproj = VectLambdaDCA.Mag();

          Minv = (PartPos1+PartNeg1).M();

          float propx = MomClosest1.x();
          float propy = MomClosest1.y();
          float propz = MomClosest1.z();
          float propE = sqrt(propx*propx + propy*propy + propz*propz + MassProton*MassProton);
          TVector3 MomProt(propx,propy,propz);
          TLorentzVector Proton4Mom(propx,propy,propz,propE);

          float pimpx = MomClosest2.x();
          float pimpy = MomClosest2.y();
          float pimpz = MomClosest2.z();
          float pimpE = sqrt(pimpx*pimpx + pimpy*pimpy + pimpz*pimpz + MassPion*MassPion);
          TVector3 MomPion(pimpx,pimpy,pimpz);
          TLorentzVector Pion4Mom(pimpx,pimpy,pimpz,pimpE);

          float pxt = propx + pimpx;
          float pyt = propy + pimpy;
          float pzt = propz + pimpz;
          float energyt = pimpE + propE;
          TVector3 MomLam(pxt,pyt,pzt);
          float LamEnergy = sqrt(LambdaDecMom.Mag2() + MassLambda*MassLambda);
          TLorentzVector Lambda4Mom(LambdaDecMom.x(),LambdaDecMom.y(),LambdaDecMom.z(),LamEnergy);

          float PV0 = sqrt(pow(pxt,2)+pow(pyt,2)+pow(pzt,2));
          LambdaMass = sqrt(energyt*energyt - pxt*pxt - pyt*pyt - pzt*pzt);
	  mLamDaughterDCA = particleDCA;
	  mLamDecLength = LambdaPrimDiff.Mag(); 
	  mLamMass = LambdaMass;
	  mLamPointAway = LambdaPrimDiff*LambdaDecMom;
	  mLamMom = MomLam.Perp();
	  mLamStart = PointLambdaToVtx;
	  if (LambdaCut(LambdaDCA, particleDCA, LambdaPrimDiff.Mag(), LambdaPrimDiff*LambdaDecMom, MomLam.Perp(), LambdaMass)) {
		if (isParticle == true) {
			StFemtoV0 vZero(pionId, protonId, isParticle, LambdaPrimDiff.Mag(), 
					particleDCA, LambdaDCA, lengths.first, lengths.second, MomLam, LambdaMass);
			vzero = vZero;}
		else	{
			StFemtoV0 vZero(pionId, protonId, isParticle, LambdaPrimDiff.Mag(),
                                        particleDCA, LambdaDCA, lengths.first, lengths.second, MomLam, LambdaMass);
			vzero = vZero;} 
			
          }
          return vzero;
}//V0Builder

StFemtoXi StPicoV0Finder::XiBuilder(int pionId, StFemtoV0 v0, bool isParticle) {
  StFemtoXi xi;
  float DCAtoPV1, DCAtoPV2;  //dca vector to primary vertex
  float px1, px2, px3, py1, py2, py3, pz1, pz2, pz3, P1, P2, P3, E1, E2, E3;
  std::pair<double,double> lengths;
  TVector3 closest1, closest2, LambdaDecMom;
  double particleDCA;
  TVector2 KaonPhiVec;
  float KaonPhi;
  TVector3 v0PointDecay, MomClosest1, MomClosest2, PointLambdaToVtx, mDca;
  double LambdaDCA;
  float Minv, LambdaMass, KinkMass, AntiKinkMass;

        StPicoTrack *mGlobTrack;
        mGlobTrack = dst->track(pionId); //Getting Positive Particle
        if(!mGlobTrack) return xi;
        if(!PionCut(mGlobTrack)) return xi;
        px1 = mGlobTrack->gMom().x();
        py1 = mGlobTrack->gMom().y();
        pz1 = mGlobTrack->gMom().z();
        TVector3 pvec1(px1,py1,pz1);
        TVector3 origin1(mGlobTrack->origin().x(),
                                      mGlobTrack->origin().y(),
                                      mGlobTrack->origin().z());
        P1 = mGlobTrack->gMom().Mag();
        E1 = sqrt(P1*P1 + MassPion*MassPion);
        StPicoPhysicalHelix helix1(pvec1, origin1, event->bField()/pow(10,14) , mGlobTrack->charge());
        DCAtoPV1 = mGlobTrack->gDCA(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());

        TLorentzVector PartPos1(px1, py1, pz1, E1);
	/*
          mGlobTrack = dst->track(pionId ); //Getting Negative Particle
          if(!mGlobTrack) return vzero;
          if (!PionCut(mGlobTrack)) return vzero;
          px2 = mGlobTrack->gMom().x();
          py2 = mGlobTrack->gMom().y();
          pz2 = mGlobTrack->gMom().z();
          P2 = mGlobTrack->gMom().Mag();
          E2 = sqrt(P2*P2 + MassPion*MassPion);
          TVector3 pvec2(px2,py2,pz2);
          TVector3 origin2(mGlobTrack->origin().x(),
                           mGlobTrack->origin().y(),
                           mGlobTrack->origin().z());
          StPicoPhysicalHelix helix2(pvec2,origin2,event->bField()/pow(10,14), mGlobTrack->charge());
          DCAtoPV2 = mGlobTrack->gDCA(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
          TLorentzVector PartNeg1(px2, py2, pz2, E2);
	*/
	  TVector3 primvertex(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
	  TLorentzVector PartNeg1 = v0.FourMom(MassLambda);
	  TVector3 LambdaOrigin = mLamStart;// decayPoint(v0, primvertex);
	  StPicoPhysicalHelix helix2(PartNeg1.Vect()*pow(PartNeg1.Vect().Perp(),-1)*100, LambdaOrigin, event->bField()/pow(10,14), 1);
          TVector3 closestDiff12;
          lengths = helix1.pathLengths(helix2);         //dca of one helix to the other
          closest1 = helix1.at(lengths.first);
          closest2 = helix2.at(lengths.second);
          closestDiff12 = closest2 - closest1;
          particleDCA = closestDiff12.Mag();
	  if (!XiDaughterCut(mLamDca, mLamDaughterDCA, mLamDecLength, mLamPointAway, mLamMom, mLamMass)) return xi;

          v0PointDecay = (closest1 + closest2)*0.5; // Calculating DCA to PV
          MomClosest1 = helix1.momentumAt(lengths.first, event->bField()/pow(10,14));
          MomClosest2 = PartNeg1.Vect();
          TVector3 XiDecMom = (MomClosest1 + MomClosest2);


//          TVector3 primvertex(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
          TVector3 LambdaPrimDiff = v0PointDecay - primvertex;
	  int XiCharge;
	  if (isParticle == true) XiCharge = -1;
	  else XiCharge = 1;
	  StPicoPhysicalHelix XiHelix(-XiDecMom, v0PointDecay, event->bField()/pow(10,14), XiCharge);
	  double DecayLength = XiHelix.pathLength(primvertex);
          TVector3 XiDCA = XiHelix.at(DecayLength);
          TVector3 XiMomAtDCA = XiHelix.momentumAt(DecayLength, event->bField()/pow(10,14));
	  TVector3 PrimDistance = v0PointDecay-primvertex;

          Minv = (PartPos1+PartNeg1).M();

          float propx = MomClosest1.x();
          float propy = MomClosest1.y();
          float propz = MomClosest1.z();
          float propE = sqrt(propx*propx + propy*propy + propz*propz + MassPion*MassPion);
          TVector3 MomProt(propx,propy,propz);
          TLorentzVector Proton4Mom(propx,propy,propz,propE);

          float pimpx = MomClosest2.x();
          float pimpy = MomClosest2.y();
          float pimpz = MomClosest2.z();
          float pimpE = sqrt(pimpx*pimpx + pimpy*pimpy + pimpz*pimpz + MassLambda*MassLambda);
          TVector3 MomPion(pimpx,pimpy,pimpz);
          TLorentzVector Pion4Mom(pimpx,pimpy,pimpz,pimpE);

          float pxt = propx + pimpx;
          float pyt = propy + pimpy;
          float pzt = propz + pimpz;
          float energyt = pimpE + propE;
          TVector3 MomXi(pxt,pyt,pzt);
          float XiEnergy = sqrt(XiDecMom.Mag2() + MassXi*MassXi);
          TLorentzVector Xi4Mom(XiDecMom.x(), XiDecMom.y(), XiDecMom.z(), XiEnergy);

          float PV0 = sqrt(pow(pxt,2)+pow(pyt,2)+pow(pzt,2));
          float XiMass = sqrt(energyt*energyt - pxt*pxt - pyt*pyt - pzt*pzt);
	 // cout << "Xi candidate: " << XiMass << " "<< (XiDCA-primvertex).Mag() << " " << DecayLength << " " << (v0PointDecay-XiMomAtDCA)*XiDecMom << endl;
//	cout << (v0PointDecay - primvertex)*XiDecMom/((v0PointDecay - primvertex).Mag()*XiDecMom.Mag()) << endl;
        if ( (decayPoint(v0, primvertex) - v0PointDecay)*MomPion > 0 && 
	     (v0PointDecay - primvertex)*XiDecMom > 0 && 
	     (decayPoint(v0, primvertex) - primvertex)*MomPion > 0 /*&&
	     (v0PointDecay - primvertex)*XiDecMom/((v0PointDecay - primvertex).Mag()*XiDecMom.Mag()) < 0.3*/  ) {
          if (XiCut((XiDCA-primvertex).Mag(), particleDCA, /*PrimDistance.Mag() */ DecayLength, 
		(v0PointDecay-primvertex)*XiDecMom, MomXi.Perp(), XiMass)) {
                        StFemtoXi XI(pionId, v0, isParticle, DecayLength,
                                        particleDCA, (XiDCA-primvertex).Mag(), lengths.first, lengths.second, MomXi, XiMass);
                        xi = XI;}}
          return xi;
}//V0Builder


bool StPicoV0Finder::PionCut(StPicoTrack *pion) {
/*	if(abs(pion->nSigmaPion()) > 3. ||
           pion->gDCA(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()) < 1. ||
           pion->gMom().Perp() < 0.15 || pion->gMom().Perp() > 10. || 
           abs(pion->gMom().Eta()) > 1.) return false;
        else return true;*/
	if (pion->gDCA(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()) < 1.) return false;
	else return true;
}

bool StPicoV0Finder::ProtonCut(StPicoTrack *proton) {
/*	if(abs(proton->nSigmaProton()) > 3. ||
           proton->gDCA(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()) < 0.3 ||
           proton->gMom().Perp() < 0.15 || proton->gMom().Perp() > 10. || 
           abs(proton->gMom().Eta()) > 1.) return false;
	else return true;*/
	if (proton->gDCA(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z()) < 0.3) return false;
	else return true;
}

bool StPicoV0Finder::LambdaCut(float DCAtoPV, float DaughterDCA, float DecayLength, float PointAway, float MomPerp, float LamMass) {
//	cout << DCAtoPV << " " <<  DaughterDCA << " " << DecayLength << " " << PointAway << endl;
	if (DCAtoPV > mDCAtoPVCut || DaughterDCA > mDaughterDCACut || DecayLength < mDecayLengthCut || PointAway < 0 ||
		 MomPerp < 0.15 || MomPerp > 10.) return false;
	else return true;
	return true;
}

bool StPicoV0Finder::XiDaughterCut(float DCAtoPV, float DaughterDCA, float DecayLength, float PointAway, float MomPerp, float LamMass) {
     /*   if (DCAtoPV > mXiDCAtoPVTopCut || DCAtoPV < mXiDCAtoPVLowCut || DaughterDCA > mXiDaughterDCACut || DecayLength < mXiDecayLengthCut || PointAway < 0 ||
                 MomPerp < 0.15 || MomPerp > 10. || fabs(LamMass-MassLambda) > mXiMassWindow) return false;
        else*/ return true;
}

bool StPicoV0Finder::XiCut(float DCAtoPV, float DaughterDCA, float DecayLength, float PointAway, float MomPerp, float XiMass) {
/*        if (DCAtoPV > mDCAtoPVCut || DaughterDCA > mDaughterDCACut || DecayLength < mDecayLengthCut || PointAway < 0 ||
                 MomPerp < 0.15 || MomPerp > 10. || fabs(XiMass-MassXi) > 0.1) return false;
        else*/ return true;
}

TVector3 StPicoV0Finder::decayPoint(StFemtoV0 v0, TVector3 primvertex) {
	TVector3 first = dst->track(v0.posDaughterId())->helix(event->bField()).at(v0.posPathLength());//pow(10,14)
	TVector3 second = dst->track(v0.negDaughterId())->helix(event->bField()).at(v0.negPathLength());
	TVector3 v0PointDecay = (first+second)*0.5;
/*	TVector3 backdir = v0.FourMom(MassLambda).Vect()*(1/(v0.FourMom(MassLambda).Vect().Mag()));
 //       TVector3 primvertex(event->primaryVertex().x(),event->primaryVertex().y(),event->primaryVertex().z());
        TVector3 LambdaPrimDiff = v0PointDecay - primvertex;

          double MinFac = - ((LambdaPrimDiff.x())*(backdir.x())
                         + (LambdaPrimDiff.y())*(backdir.y())
                         + (LambdaPrimDiff.z())*(backdir.z()))/(backdir.Mag());
          TVector3 PointLambdaToVtx = v0PointDecay + MinFac*backdir;
	  return PointLambdaToVtx;*/
	return v0PointDecay;
}

StPicoPhysicalHelix StPicoV0Finder::helix(StFemtoV0 v0, float field) {
	StPicoPhysicalHelix V0Helix(v0.FourMom(MassLambda).Vect()*pow(v0.FourMom(MassLambda).Vect().Perp(),-1)*100, decayPoint(v0,TVector3(0.,0.,0.)), event->bField()/pow(10,14), 1.);
	return V0Helix;	
}

