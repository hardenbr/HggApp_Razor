#include "cajun/json/reader.h"
#include "cajun/json/elements.h"

#include <fstream>
#include <sstream>
#include <math.h>

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "EgammaAnalysisTools/include/Egamma.h"

using namespace bits;
using namespace std;

Egamma::Egamma(TTree *tree) : EgammaBase(tree)
{
  jsonFile = "";
  lastFile = "";

  // H->ZZ/WW effective areas computed in 2012
  float Aeff_neu_dr04_[7] = { 0.045, 0.065, 0.068, 0.057, 0.058, 0.061, 0.110 };
  float Aeff_pho_dr04_[7] = { 0.140, 0.130, 0.079, 0.130, 0.150, 0.160, 0.180 };
  
  // H->ZZ 2011 detector based effective areas
  float Aeff_ecal_dr03_[2] = { 0.078, 0.046 };
  float Aeff_hcal_dr03_[2] = { 0.026, 0.072 };
  
  for(int i=0;i<7;i++) {
    Aeff_neu_dr04[i]= Aeff_neu_dr04_[i];
    Aeff_pho_dr04[i]= Aeff_pho_dr04_[i];
  }
  for(int i=0;i<2;i++) {
    Aeff_ecal_dr03[i]= Aeff_ecal_dr03_[i];
    Aeff_hcal_dr03[i]= Aeff_hcal_dr03_[i];
  }

}

Egamma::~Egamma()
{
  // By this time, the destructor of EgammaBase has not yet been called.
  // This means that the tree has not yet been deleted.
  // So, we do nothing here.
}

void Egamma::setRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  requiredTriggers=reqTriggers;
}

bool Egamma::hasPassedHLT() {
  Utils anaUtils;
  return anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
}

void Egamma::setJsonGoodRunList(const std::string& jsonFilePath)
{
  jsonFile=jsonFilePath;
}

void Egamma::fillRunLSMap()
{
  
  if (jsonFile == "")
    {
      std::cout << "Cannot fill RunLSMap. json file not configured" << std::endl;
      return;
    }

  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile.c_str());
  if (!jsonFileStream.is_open())
    {
      std::cout << "Unable to open file " << jsonFile << std::endl;
      return;
    }

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun)
    {
      const json::Array& lsSegment = (*itRun).element;
      LSSegments thisRunSegments; 
      for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator)
	{
	  json::Array lsSegment=(*lsIterator);
	  json::Number lsStart=lsSegment[0];	   
	  json::Number lsEnd=lsSegment[1];
	  aLSSegment thisSegment;
	  thisSegment.first=lsStart.Value();
	  thisSegment.second=lsEnd.Value();
	  thisRunSegments.push_back(thisSegment);
	  //	   std::pair<int, int> lsSegment=std::pair<int, int>(atoi(,lsIterator[1]); 
	}
      goodRunLS.insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
    }


  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS.size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS.begin(); itR!=goodRunLS.end(); ++itR)
    {
      std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
      for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
	std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] "; 
      std::cout << std::endl;
    }
}

bool Egamma::isGoodRunLS()
{
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);
  if (thisRun == goodRunLS.end())
    return false;
  //  std::cout << runNumber << " found in the good run map" << std::endl;
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg)
    {
      //      std::cout << "Range is [" << (*iSeg).first << "," << (*iSeg).second << "]" << std::endl;
      if ( lumiBlock >= (*iSeg).first && lumiBlock <= (*iSeg).second)
	return true;
    }
  return false;
}

bool Egamma::reloadTriggerMask(bool newVersion)
{
  if(newVersion) {
    std::vector<int> triggerMask;
    for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
      {   
        for(unsigned int i=0; i<nameHLT->size(); i++)
          {
            // if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
            if(nameHLT->at(i).find(*fIter) != std::string::npos)
              {
                triggerMask.push_back( indexHLT[i] ) ;
                break;
              }
          }
      }
    m_requiredTriggers = triggerMask;
  } else {
    TString fileName=((TChain*)fChain)->GetFile()->GetName();
    if ( TString(lastFile) != fileName )
      {

        std::cout << "[ReloadTriggerMask]::File has changed reloading trigger mask" << std::endl;
        lastFile = fileName;
        TTree *treeCond;
        std::cout << "[ReloadTriggerMask]::Opening " << fileName << std::endl;
        treeCond = (TTree*)((TChain*)fChain)->GetFile()->Get("Conditions");
        int           nHLT_;
        std::vector<std::string>  *nameHLT_;
        std::vector<unsigned int> *indexHLT_;

        //To get the pointers for the vectors
        nameHLT_=0;
        indexHLT_=0;

        treeCond->SetBranchAddress("nHLT", &nHLT_);
        treeCond->SetBranchAddress("nameHLT", &nameHLT_);
        treeCond->SetBranchAddress("indexHLT", &indexHLT_);
        treeCond->GetEntry(0);

        std::vector<int> triggerMask;
        for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
          {
            for(unsigned int i=0; i<nameHLT_->size(); i++) 
              {
                if( !strcmp ((*fIter).c_str(), nameHLT_->at(i).c_str() ) ) 
                  {
                    triggerMask.push_back( indexHLT_->at(i) ) ;
                    break;
                  }
              }
          }
        m_requiredTriggers = triggerMask;
        for (int i=0;i<m_requiredTriggers.size();++i)
          std::cout << "[ReloadTriggerMask]::Requiring bit " << m_requiredTriggers[i] << " " << requiredTriggers[i] << std::endl;
      }
  }
}

float Egamma::mT3(TLorentzVector pl1, TLorentzVector pl2, TVector3 met) {
  float pTll = (pl1.Vect() + pl2.Vect()).Pt();
  float mll = (pl1 + pl2).M();
  float El = sqrt(pTll*pTll + mll*mll);
  float pTnu = met.Pt();
  float Enu = sqrt(pTnu*pTnu + mll*mll);
  float Ex = (pl1+pl2).X() + met.X();
  float Ey = (pl1+pl2).Y() + met.Y();
  float mnu = mll;

  return sqrt(mll*mll + mnu*mnu + 2*(El*Enu-Ex*Ex-Ey*Ey));
}

float Egamma::likelihoodRatio(int eleIndex, ElectronLikelihood &lh) {
  LikelihoodMeasurements measurements;
  Utils anaUtils;
  bool inEB=anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex],isEB);
  measurements.pt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);
  if(inEB && fabs(etaEle[eleIndex])<1.0) measurements.subdet = 0;
  else if (inEB && fabs(etaEle[eleIndex])>=1.0) measurements.subdet = 1;
  else measurements.subdet = 2;
  measurements.deltaPhi = deltaPhiAtVtxEle[eleIndex];
  measurements.deltaEta = deltaEtaAtVtxEle[eleIndex];
  measurements.eSeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  measurements.eSuperClusterOverP = eSuperClusterOverPEle[eleIndex];
  measurements.hadronicOverEm = hOverEEle[eleIndex];
  measurements.sigmaIEtaIEta = SigmaiEiE(eleIndex);
  measurements.sigmaIPhiIPhi = SigmaiPiP(eleIndex);
  measurements.fBrem = fbremEle[eleIndex];
  measurements.nBremClusters = nbremsEle[eleIndex];
  int gsftrack = gsfTrackIndexEle[eleIndex];
  TVector3 pIn(pxGsfTrack[gsftrack],pyGsfTrack[gsftrack],pzGsfTrack[gsftrack]);
  measurements.OneOverEMinusOneOverP = 1./(eSuperClusterOverPEle[eleIndex]*pIn.Mag()) - 1./pIn.Mag();
  return lh.resultLog(measurements);
}

/// sigma ieta ieta of the seed cluster (ecal-driven/tracker-driven)
float Egamma::SigmaiEiE(int electron) {
  float see;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    see = sqrt(covIEtaIEtaSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      see = sqrt(covIEtaIEtaPFSC[sc]);
    } else {
      see = 999.;
    }
  }
  return see;
}

/// sigma iphi iphi of the seed cluster (ecal-driven/tracker-driven)
float Egamma::SigmaiPiP(int electron) {
  float spp;
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[electron], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[electron];
    spp = sqrt(covIPhiIPhiSC[sc]);
  } else {
    int sc = PFsuperClusterIndexEle[electron];
    if(sc>-1) {
      spp = sqrt(covIPhiIPhiPFSC[sc]);
    } else {
      spp = 999.;
    }
  }
  return spp;
}

bool Egamma::triggerMatch(float eta, float phi, float Dr){

  bool match=false;
  for( int i=0; i<m_requiredTriggers.size(); i++ ) {  // loop over require trigger paths
    
    int pathIndex=m_requiredTriggers[i];
    // std::cout << "testing trigger " << pathIndex << " with " << sizePassing[pathIndex] << " passing objects" << std::endl; 
    
    if( sizePassing[pathIndex]>  0 ) {  //some object has passed the required trigger 
      
      for(int np = 0; np < sizePassing[pathIndex]; np++ ){
        int iP = indexPassing[ indexPassingPerPath[pathIndex] +np];
        // std::cout << "passing object eta: " << triggerObsEta[iP] << " phi: " <<  triggerObsPhi[iP] << std::endl; 

        if(DeltaR(eta, phi,triggerObsEta[iP],  triggerObsPhi[iP] ) < Dr){
          match=true;
          //std::cout << "MATCH!" <<std::endl;	
          break;
        }
      }
    }
    if(match)  //it's enough if one path matches	
      break;
  }
  return match;
}

bool Egamma::triggerMatchThreshold(float pt){

  bool match=false;

  for( int i=0; i<m_requiredTriggers.size(); i++ ) {  // loop over require trigger paths
    
    int pathIndex=m_requiredTriggers[i];
    if( sizePassing[pathIndex]>  0 ) {  //some object has passed the required trigger 
      
      for(int np = 0; np < sizePassing[pathIndex]; np++ ){
        int iP = indexPassing[ indexPassingPerPath[pathIndex] +np];

	if (triggerObsPt[iP]>pt) {
          match=true;
          break;
        }
      }
    }
    if(match)  //it's enough if one path matches	
      break;
  }
  return match;
}

// two highest pT electrons 
std::pair<int,int> Egamma::getBestGoodElePair(std::vector<int> goodElectrons) {
  
  int theEle1=-1;
  int theEle2=-1;
  float maxPt1=-1000.;
  float maxPt2=-1001.;
  for(int iEle=0;iEle<goodElectrons.size();iEle++) {
    int eleIndex = goodElectrons[iEle];
    TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
    float thisPt=pEle.Pt();
    if (thisPt>maxPt1 && thisPt>maxPt2){ maxPt2 = maxPt1; maxPt1 = thisPt; theEle2 = theEle1; theEle1 = eleIndex; }
    if (thisPt<maxPt1 && thisPt>maxPt2){ maxPt2 = thisPt; theEle2 = eleIndex; }
  }
  return std::make_pair(theEle1,theEle2);
}

// dxy parameter with respect to PV for tracks
double Egamma::trackDxyPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  return ( - (trackVPos.X()-PVPos.X())*trackMom.Y() + (trackVPos.Y()-PVPos.Y())*trackMom.X() ) / trackMom.Pt(); 
}

/// dz parameter with respect to PV for tracks
double Egamma::trackDzPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  float trackPt = trackMom.Pt();
  return (trackVPos.Z()-PVPos.Z()) - ((trackVPos.X()-PVPos.X())*trackMom.X()+(trackVPos.Y()-PVPos.Y())*trackMom.Y())/trackPt *trackMom.Pz()/trackPt; 
}

/// dsz parameter with respect to PV for tracks
double Egamma::trackDszPV(TVector3 PVPos, TVector3 trackVPos, TVector3 trackMom) {
  float trackPt = trackMom.Pt();
  float trackP  = trackMom.Mag();
  return (trackVPos.Z()-PVPos.Z())*trackPt/trackP - ((trackVPos.X()-PVPos.X())*trackMom.X()+(trackVPos.Y()-PVPos.Y())*trackMom.Y())/trackPt *trackMom.Pz()/trackP; 
}

/// dxy, dz and dsz parameters with respect to PV for electrons
double Egamma::eleDxyPV(int iele, int iPV) {
  TVector3 PVPos(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  int gsfTrack = gsfTrackIndexEle[iele];
  TVector3 lepVPos(trackVxGsfTrack[gsfTrack],trackVyGsfTrack[gsfTrack],trackVzGsfTrack[gsfTrack]);
  TVector3 lepMom(pxEle[iele],pyEle[iele],pzEle[iele]);
  return trackDxyPV(PVPos,lepVPos,lepMom);
}

double Egamma::eleDzPV(int iele, int iPV) {
  TVector3 PVPos(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  int gsfTrack = gsfTrackIndexEle[iele];
  TVector3 lepVPos(trackVxGsfTrack[gsfTrack],trackVyGsfTrack[gsfTrack],trackVzGsfTrack[gsfTrack]);
  TVector3 lepMom(pxEle[iele],pyEle[iele],pzEle[iele]);
  return trackDzPV(PVPos,lepVPos,lepMom);
}

double Egamma::eleDszPV(int iele, int iPV) {
  TVector3 PVPos(PVxPV[iPV],PVyPV[iPV],PVzPV[iPV]);
  int gsfTrack = gsfTrackIndexEle[iele];
  TVector3 lepVPos(trackVxGsfTrack[gsfTrack],trackVyGsfTrack[gsfTrack],trackVzGsfTrack[gsfTrack]);
  TVector3 lepMom(pxEle[iele],pyEle[iele],pzEle[iele]);
  return trackDszPV(PVPos,lepVPos,lepMom);
}

// using HWW BDT (paper 2011) 
float Egamma::eleBDT(ElectronIDMVA *mva, int eleIndex) {
  
  if(mva==0) {
    std::cout << "electron BDT not created/initialized. BIG PROBLEM. Returning false output -999!" << std::endl; 
    return -999.;
  }
  
  int gsfTrack = gsfTrackIndexEle[eleIndex]; 
  double gsfsign   = (-eleDxyPV(eleIndex,0) >=0 ) ? 1. : -1.;

  float ElePt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);

  float EleDEtaIn = deltaEtaAtVtxEle[eleIndex];
  float EleDPhiIn = deltaPhiAtVtxEle[eleIndex];
  float EleHoverE = hOverEEle[eleIndex];
  float EleD0 = gsfsign * transvImpactParGsfTrack[gsfTrack];
  float EleFBrem = fbremEle[eleIndex];
  float EleEOverP = eSuperClusterOverPEle[eleIndex];
  float EleESeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  float EleNBrem = nbremsEle[eleIndex];
  TVector3 pInGsf(pxGsfTrack[gsfTrack],pyGsfTrack[gsfTrack],pzGsfTrack[gsfTrack]);

  float EleIP3d = gsfsign * impactPar3DGsfTrack[gsfTrack];
  float EleIP3dSig = EleIP3d/impactPar3DErrorGsfTrack[gsfTrack];

  // we have not pout and seed cluster energy in the trees. Gymnastyc...
  float Pout = pInGsf.Mag() - fbremEle[eleIndex] * pInGsf.Mag();
  float SCSeedEnergy = EleESeedClusterOverPout * Pout;
  float EleESeedClusterOverPIn = SCSeedEnergy/pInGsf.Mag();      

  float EleSigmaIEtaIEta, EleSigmaIPhiIPhi, EleOneOverEMinusOneOverP, EleSCEta;

  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    EleSigmaIEtaIEta = sqrt(covIEtaIEtaSC[sc]);
    EleSigmaIPhiIPhi = sqrt(covIPhiIPhiSC[sc]);
    EleOneOverEMinusOneOverP = 1./energySC[sc]  - 1./pInGsf.Mag();
    EleSCEta = etaSC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      EleSigmaIEtaIEta = sqrt(covIEtaIEtaPFSC[sc]);
      EleSigmaIPhiIPhi = sqrt(covIPhiIPhiPFSC[sc]);
      EleOneOverEMinusOneOverP = 1./energyPFSC[sc]  - 1./pInGsf.Mag();
      EleSCEta = etaPFSC[sc];
    } else {
      EleSigmaIEtaIEta = 999.;
      EleSigmaIPhiIPhi = 999.;
      EleOneOverEMinusOneOverP = 999.;
      EleESeedClusterOverPIn = 999.;
      EleSCEta = 0.;
    }
  }

  return mva->MVAValue(ElePt, EleSCEta,
                       EleSigmaIEtaIEta,
                       EleDEtaIn,
                       EleDPhiIn,
                       EleHoverE,
                       EleD0,
                       EleFBrem,
                       EleEOverP,
                       EleESeedClusterOverPout,
                       EleSigmaIPhiIPhi,
                       EleNBrem,
                       EleOneOverEMinusOneOverP,
                       EleESeedClusterOverPIn,
                       EleIP3d,
                       EleIP3dSig );

}

// Si proposal for 2012
float Egamma::eleBDTWithIso(ElectronIDMVA *mva, int eleIndex) {
  
  if(mva==0) {
    std::cout << "electron BDT not created/initialized. BIG PROBLEM. Returning false output -999!" << std::endl; 
    return -999.;
  }
  
  int gsfTrack = gsfTrackIndexEle[eleIndex]; 
  double gsfsign   = (-eleDxyPV(eleIndex,0) >=0 ) ? 1. : -1.;

  float ElePt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);

  float EleDEtaIn = deltaEtaAtVtxEle[eleIndex];
  float EleDPhiIn = deltaPhiAtVtxEle[eleIndex];
  float EleDEtaCalo = deltaEtaAtCaloEle[eleIndex];
  float EleDPhiCalo = deltaPhiAtCaloEle[eleIndex];
  float EleHoverE = hOverEEle[eleIndex];
  float EleD0 = gsfsign * transvImpactParGsfTrack[gsfTrack];
  float EleFBrem = fbremEle[eleIndex];
  float EleEOverP = eSuperClusterOverPEle[eleIndex];
  float EleESeedClusterOverPout = eSeedOverPoutEle[eleIndex];
  float EleNBrem = nbremsEle[eleIndex];
  float EleGSFChi2 = trackNormalizedChi2GsfTrack[gsfTrack];

  TVector3 pInGsf(pxGsfTrack[gsfTrack],pyGsfTrack[gsfTrack],pzGsfTrack[gsfTrack]);

  float EleIP3d = gsfsign * impactPar3DGsfTrack[gsfTrack];
  float EleIP3dSig = EleIP3d/impactPar3DErrorGsfTrack[gsfTrack];

  // we have not pout and seed cluster energy in the trees. Gymnastyc...
  float Pout = pInGsf.Mag() - fbremEle[eleIndex] * pInGsf.Mag();
  float SCSeedEnergy = EleESeedClusterOverPout * Pout;
  float EleESeedClusterOverPIn = SCSeedEnergy/pInGsf.Mag();      

  float EleSigmaIEtaIEta, EleSigmaIPhiIPhi, EleOneOverEMinusOneOverP, EleSCEta, EleR9, 
    EleCovIEtaIPhi, EleEtaWidth, ElePhiWidth, EleESoRaw;

  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    EleSigmaIEtaIEta = sqrt(covIEtaIEtaSC[sc]);
    EleSigmaIPhiIPhi = sqrt(covIPhiIPhiSC[sc]);
    EleCovIEtaIPhi = covIEtaIPhiSC[sc];
    EleOneOverEMinusOneOverP = 1./energySC[sc]  - 1./pInGsf.Mag();
    EleEtaWidth = etaWidthSC[sc];
    ElePhiWidth = phiWidthSC[sc];
    EleR9 = e3x3SC[sc]/rawEnergySC[sc];
    EleESoRaw = esEnergySC[sc]/rawEnergySC[sc];
    EleSCEta = etaSC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      EleSigmaIEtaIEta = sqrt(covIEtaIEtaPFSC[sc]);
      EleSigmaIPhiIPhi = sqrt(covIPhiIPhiPFSC[sc]);
      EleCovIEtaIPhi = covIEtaIPhiPFSC[sc];
      EleOneOverEMinusOneOverP = 1./energyPFSC[sc]  - 1./pInGsf.Mag();
      EleEtaWidth = etaWidthPFSC[sc];
      ElePhiWidth = phiWidthPFSC[sc];
      EleR9 = e3x3PFSC[sc]/rawEnergyPFSC[sc];
      EleESoRaw = esEnergyPFSC[sc]/rawEnergyPFSC[sc];
      EleSCEta = etaPFSC[sc];
    } else {
      EleSigmaIEtaIEta = 999.;
      EleSigmaIPhiIPhi = 999.;
      EleOneOverEMinusOneOverP = 999.;
      EleESeedClusterOverPIn = 999.;
      EleEtaWidth = 999.;
      ElePhiWidth = 999.;
      EleESoRaw = esEnergyPFSC[sc]/rawEnergyPFSC[sc];
      EleR9 = 999.;
      EleSCEta = 0.;
    }
  }

  float EleChargedIso03OverPt = (pfCandChargedIso03Ele[eleIndex] - rhoFastjet * SiElectronEffectiveArea(kEleChargedIso03,EleSCEta))/ElePt;
  float EleGammaIso03OverPt = (pfCandPhotonIso03Ele[eleIndex] 
			       - rhoFastjet * SiElectronEffectiveArea(kEleGammaIso03,EleSCEta)
			       + rhoFastjet * SiElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,EleSCEta))/ElePt;
  float EleNeutralHadronIso03OverPt = (pfCandNeutralIso03Ele[eleIndex] 
				       - rhoFastjet * SiElectronEffectiveArea(kEleNeutralHadronIso03,EleSCEta)
				       + rhoFastjet * SiElectronEffectiveArea(kEleNeutralHadronIso007,EleSCEta))/ElePt;

  float EleChargedIso04OverPt = (pfCandChargedIso04Ele[eleIndex] - rhoFastjet * SiElectronEffectiveArea(kEleChargedIso04,EleSCEta))/ElePt;
  float EleGammaIso04OverPt = (pfCandPhotonIso04Ele[eleIndex] 
			       - rhoFastjet * SiElectronEffectiveArea(kEleGammaIso04,EleSCEta)
			       + rhoFastjet * SiElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,EleSCEta))/ElePt;
  float EleNeutralHadronIso04OverPt = (pfCandNeutralIso04Ele[eleIndex] 
				       - rhoFastjet * SiElectronEffectiveArea(kEleNeutralHadronIso04,EleSCEta)
				       + rhoFastjet * SiElectronEffectiveArea(kEleNeutralHadronIso007,EleSCEta))/ElePt;

  return mva->MVAValueWithIso(ElePt, EleSCEta,
			      EleSigmaIEtaIEta,
			      EleDEtaIn,
			      EleDPhiIn,
			      EleHoverE,
			      EleD0,
			      EleFBrem,
			      EleEOverP,
			      EleESeedClusterOverPout,
			      EleSigmaIPhiIPhi,
			      EleNBrem,
			      EleOneOverEMinusOneOverP,
			      EleESeedClusterOverPIn,
			      EleIP3d,
			      EleIP3dSig,
			      EleGSFChi2,
			      EleDEtaCalo,
			      EleDPhiCalo,
			      EleR9,
			      EleEtaWidth,
			      ElePhiWidth,
			      EleCovIEtaIPhi,
			      EleESoRaw,
			      EleChargedIso03OverPt,
			      EleNeutralHadronIso03OverPt,
			      EleGammaIso03OverPt,
			      EleChargedIso04OverPt,
			      EleNeutralHadronIso04OverPt,
			      EleGammaIso04OverPt);

}

// using HZZ BDT
float Egamma::eleBDT(ElectronIDMVAHZZ *mva, int eleIndex) {
  
  if(mva==0) {
    std::cout << "electron BDT not created/initialized. BIG PROBLEM. Returning false output -999!" << std::endl; 
    return -999.;
  }
  
  TLorentzVector probeP4;
  probeP4.SetXYZT(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex],energyEle[eleIndex]);

  int gsfTrack = gsfTrackIndexEle[eleIndex]; 
  int kfTrack = trackIndexEle[eleIndex];
  float ElePt = GetPt(pxEle[eleIndex],pyEle[eleIndex]);

  float EleFBrem = fbremEle[eleIndex];
  float EleNBrems = nbremsEle[eleIndex];
  float EleDEtaIn = deltaEtaAtVtxEle[eleIndex];
  float EleDPhiIn = deltaPhiAtVtxEle[eleIndex];
  float EleDEtaCalo = deltaEtaAtCaloEle[eleIndex];
  float EleDPhiCalo = deltaPhiAtCaloEle[eleIndex];
  float EleHoverE = hOverEEle[eleIndex];
  float EleSuperClusterEOverP = eSuperClusterOverPEle[eleIndex];
  float EleEeleOverPout = eEleClusterOverPoutEle[eleIndex];
  float EleEseedOverPout = eSeedOverPoutEle[eleIndex];
  float EleDEtaEleOut = deltaEtaEleClusterTrackAtCaloEle[eleIndex];
  float EleKFChi2 = (kfTrack>-1) ? trackNormalizedChi2Track[kfTrack] : 0.0;
  float EleGSFChi2 = trackNormalizedChi2GsfTrack[gsfTrack];
  float EleKFHits = (kfTrack>-1) ? trackerLayersWithMeasurementTrack[kfTrack] : -1.0;
  float EleMissHits = expInnerLayersGsfTrack[gsfTrack];
  float EleDistConv = convDistEle[eleIndex];
  float EleDcotConv = convDcotEle[eleIndex];
  float NVtx = float(nPV);
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  float EleEcalSeeded = (ecaldriven) ? 1. : 0.;
  float EleMatchConv = (hasMatchedConversionEle[eleIndex]) ? 1. : 0.;

  // IP variables
  double gsfsign   = (-eleDxyPV(eleIndex,0) >=0 ) ? 1. : -1.;
  float EleD0 = gsfsign * transvImpactParGsfTrack[gsfTrack];
  float EleIP3d = gsfsign * impactPar3DGsfTrack[gsfTrack];
  float EleIP3dSig = EleIP3d/impactPar3DErrorGsfTrack[gsfTrack];

  float EleSigmaIEtaIEta, EleSigmaIPhiIPhi, EleSigmaIEtaIPhi, EleE1x5E5x5, EleSCEta, EleEtaWidth, ElePhiWidth,
    EleR9, Ele1oEm1oP, EleESoRaw;

  // used only in the case of triggering electrons
  float EleHWWPresel = (float)isDenomFake_smurfs(eleIndex);

  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    EleSigmaIEtaIEta = sqrt(covIEtaIEtaSC[sc]);
    EleSigmaIPhiIPhi = sqrt(covIPhiIPhiSC[sc]);
    EleSigmaIEtaIPhi = covIEtaIPhiSC[sc]/(sqrt(covIEtaIEtaSC[sc])*sqrt(covIPhiIPhiSC[sc]));
    EleE1x5E5x5 = (e5x5SC[sc] - e1x5SC[sc])/e5x5SC[sc];
    EleEtaWidth = etaWidthSC[sc];
    ElePhiWidth = phiWidthSC[sc];
    EleR9 = e3x3SC[sc]/rawEnergySC[sc];
    Ele1oEm1oP = 1./energySC[sc]  - 1./probeP4.Vect().Mag();
    EleESoRaw = esEnergySC[sc]/rawEnergySC[sc];
    EleSCEta = etaSC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      EleSigmaIEtaIEta = sqrt(covIEtaIEtaPFSC[sc]);
      EleSigmaIPhiIPhi = sqrt(covIPhiIPhiPFSC[sc]);
      EleSigmaIEtaIPhi = covIEtaIPhiPFSC[sc]/(sqrt(covIEtaIEtaPFSC[sc])*sqrt(covIPhiIPhiPFSC[sc]));
      EleE1x5E5x5 = (e5x5PFSC[sc] - e1x5PFSC[sc])/e5x5PFSC[sc];
      EleEtaWidth = etaWidthPFSC[sc];
      ElePhiWidth = phiWidthPFSC[sc];
      EleR9 = e3x3PFSC[sc]/rawEnergyPFSC[sc];
      Ele1oEm1oP = 1./energyPFSC[sc]  - 1./probeP4.Vect().Mag();
      EleESoRaw = esEnergyPFSC[sc]/rawEnergyPFSC[sc];
      EleSCEta = etaPFSC[sc];
    } else {
      EleSigmaIEtaIEta = 999.;
      EleSigmaIPhiIPhi = 999.;
      EleSigmaIEtaIPhi = 999.;
      EleE1x5E5x5 = 999.;
      EleEtaWidth = 999.;
      ElePhiWidth = 999.;  
      EleR9 = 999.;
      Ele1oEm1oP = 999.;
      EleESoRaw = esEnergyPFSC[sc]/rawEnergyPFSC[sc];
      EleSCEta = 0.;
    }
  }

  return mva->MVAValue(ElePt, EleSCEta,
                       EleFBrem,
                       EleNBrems,
                       EleDEtaIn,
                       EleDPhiIn,
                       EleDPhiCalo,
                       EleDEtaCalo,
                       EleDEtaEleOut,
                       EleSigmaIEtaIEta,
                       EleSigmaIPhiIPhi,
                       EleSigmaIEtaIPhi,
                       EleHoverE,
                       EleSuperClusterEOverP,
                       EleE1x5E5x5,
                       EleR9,
                       EleESoRaw,
                       EleEseedOverPout,
                       EleEeleOverPout,
                       Ele1oEm1oP,
                       EleKFChi2,
                       EleKFHits,
                       EleGSFChi2,
                       EleMissHits,
                       EleDistConv,
                       EleDcotConv,
                       EleMatchConv,
                       NVtx,
                       EleEcalSeeded,
                       EleEtaWidth,
                       ElePhiWidth,
                       EleD0,
                       EleIP3d,
                       EleIP3dSig,
                       EleHWWPresel);

}

bool Egamma::passEleBDT(float pt, float eta, float bdt) {

  if(pt < 20 && fabs(eta) < 1.0) return (bdt > 0.139);
  if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdt > 0.525);
  if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdt > 0.543);
  if(pt >= 20 && fabs(eta) < 1.0) return (bdt > 0.947);
  if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdt > 0.950);
  if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdt > 0.884);

  // here we are cutting the events with |SC eta|>2.5. If the acceptance is done with |ele eta|<2.5 then will cut some event more. Fine. Synch with this.
  return false;

}

float Egamma::combPFIsoEACorr(int iele) {
  float pt = GetPt(pxEle[iele],pyEle[iele]);
  Utils anaUtils;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[iele], isEcalDriven);
  float eta = (ecaldriven) ? etaSC[superClusterIndexEle[iele]] : etaPFSC[PFsuperClusterIndexEle[iele]];

  float combIso = pfCandChargedIso04Ele[iele] + pfCandNeutralIso04Ele[iele] + pfCandPhotonIso04Ele[iele];

  // EA correction is calculated for NH+PHO only (charged Had is computed with PFnoPU): cha + (neu - EA*rho)
  if(pt <  20 && fabs(eta) <  1.479) return combIso - 0.164 * rhoFastjet;
  if(pt >= 20 && fabs(eta) <  1.479) return combIso - 0.175 * rhoFastjet;
  if(pt <  20 && fabs(eta) >= 1.479) return combIso - 0.210 * rhoFastjet;
  if(pt >= 20 && fabs(eta) >= 1.479) return combIso - 0.333 * rhoFastjet;
  return 0.0;
}

// denominator for fake rate: for HtoWW, egamma triggers
int Egamma::isDenomFake_HwwEgamma(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);

  // acceptance for the fake electron
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )        isGoodDenom = false;

  // barrel or endcap
  bool isEleEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[theEle], isEB);    

  // taking shower shape                                                                                                             
  int sc;
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  float thisSigmaIeIe = -1.;
  if (ecalDriven) {
    sc = superClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaSC[sc]);
  }
  if (!ecalDriven) {
    sc = PFsuperClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaPFSC[sc]);
  }
  if ( sc<0 ) { isGoodDenom = false; }

  // sigmaIetaIeta                                                                                                                   
  if ( isEleEB && thisSigmaIeIe>0.01) { isGoodDenom = false; }
  if (!isEleEB && thisSigmaIeIe>0.03) { isGoodDenom = false; }

  // H/E
  if ( isEleEB && hOverEEle[theEle]>0.12) isGoodDenom = false;   
  if (!isEleEB && hOverEEle[theEle]>0.10) isGoodDenom = false;
  
  // deltaEta
  if ( isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.007) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.009) ) isGoodDenom = false;
  
  // deltaPhi
  if ( isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.15) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.10) ) isGoodDenom = false;

  // isolation 
  float ecalIsol    = (dr03EcalRecHitSumEtEle[theEle])/p3Ele.Pt();
  float hcalIsol    = (dr03HcalTowerSumEtEle[theEle])/p3Ele.Pt();
  float trackerIsol = (dr03TkSumPtEle[theEle])/p3Ele.Pt();                
  if(ecalIsol>0.2)    isGoodDenom = false;
  if(hcalIsol>0.2)    isGoodDenom = false;
  if(trackerIsol>0.2) isGoodDenom = false;                                

  if(isGoodDenom) return 1;
  else return 0;
}

// denominator for fake rate: for HtoWW, egamma triggers, same as smurfs
int Egamma::isDenomFake_smurfs(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // acceptance for the fake electron
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )        isGoodDenom = false;

  // taking shower shape                                                                                                             
  int sc;
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  float thisSigmaIeIe = -1.;
  float scEta = -1.;                                                               
  if (ecalDriven) {
    sc = superClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaSC[sc]);
    scEta = etaSC[sc];
  }
  if (!ecalDriven) {
    sc = PFsuperClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaPFSC[sc]);
    scEta = etaPFSC[sc];
  }
  if ( sc<0 ) { isGoodDenom = false; }

  // barrel or endcap
  bool isEleEB = false;
  if (fabs(scEta)<1.479) isEleEB = true;   
  
  // sigmaIetaIeta                                                                                                                   
  if ( isEleEB && thisSigmaIeIe>0.01) { isGoodDenom = false; }
  if (!isEleEB && thisSigmaIeIe>0.03) { isGoodDenom = false; }

  // isolation
  float ecalIsolAbs = 0.0;
  if ( isEleEB ) ecalIsolAbs = max(0.0,dr03EcalRecHitSumEtEle[theEle]-1.0);
  else ecalIsolAbs = dr03EcalRecHitSumEtEle[theEle];
  float ecalIsol = ecalIsolAbs/p3Ele.Pt(); 
  float hcalIsol    = (dr03HcalTowerSumEtEle[theEle])/p3Ele.Pt();
  float trackerIsol = (dr03TkSumPtEle[theEle])/p3Ele.Pt();                
  if(ecalIsol>0.2)    isGoodDenom = false;
  if(hcalIsol>0.2)    isGoodDenom = false;
  if(trackerIsol>0.2) isGoodDenom = false;                                

  // H/E
  if ( isEleEB && hOverEEle[theEle]>0.12) isGoodDenom = false;   
  if (!isEleEB && hOverEEle[theEle]>0.10) isGoodDenom = false;
  
  // deltaEta
  if ( isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.007) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.009) ) isGoodDenom = false;
  
  // deltaPhi
  if ( isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.15) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.10) ) isGoodDenom = false;
  
  // full conversion rejection 
  int gsf = gsfTrackIndexEle[theEle];
  int missHits = expInnerLayersGsfTrack[gsf];
  bool matchConv = hasMatchedConversionEle[theEle];
  if (missHits>0 || matchConv) isGoodDenom = false;

  // impact parameter cuts 
  float dxyEle = transvImpactParGsfTrack[gsf];
  float dzEle  = eleDzPV(theEle,0);
  if (fabs(dxyEle)>0.02) isGoodDenom = false;
  if (fabs(dzEle)>0.10)  isGoodDenom = false;

  if(isGoodDenom) return 1;
  else return 0;
}

double Egamma::SiElectronEffectiveArea(ElectronEffectiveAreaType type, double Eta) {

  double EffectiveArea = 0;
  
  if (fabs(Eta) < 1.0) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.017;
    if (type == kEleGammaIso03) EffectiveArea = 0.045;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.014;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.034;
    if (type == kEleGammaIso04) EffectiveArea = 0.079;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.014;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00016;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00016;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.025;
    if (type == kEleGammaIso03) EffectiveArea = 0.052;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.030;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.050;
    if (type == kEleGammaIso04) EffectiveArea = 0.073;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.030;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00022;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00022;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.030;
    if (type == kEleGammaIso03) EffectiveArea = 0.170;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.134;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.060;
    if (type == kEleGammaIso04) EffectiveArea = 0.187;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.134;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00030;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00026;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00002;        
  } else if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.25 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.022;
    if (type == kEleGammaIso03) EffectiveArea = 0.623;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.516;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.055;
    if (type == kEleGammaIso04) EffectiveArea = 0.659;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.517;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00054;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00045;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00003;
  } else if (fabs(Eta) >= 2.25 && fabs(Eta) < 2.5 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.018;
    if (type == kEleGammaIso03) EffectiveArea = 1.198;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 1.049;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.073;
    if (type == kEleGammaIso04) EffectiveArea = 1.258;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 1.051;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00082;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00066;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00004;
  }
  return EffectiveArea;  
}

int Egamma::indexSeedBC(int sc, int ecaldriven) {
  std::vector<int> BCs;
  if(ecaldriven) {
    for(int ibc=0;ibc<nBC;++ibc) {
      if(indexSCBC[ibc]==sc) BCs.push_back(ibc);
    }
    // find the seed BC
    float maxE=-1.;
    int seed=-1;
    for(unsigned int ibc=0; ibc<BCs.size(); ++ibc) {
      if(energyBC[BCs[ibc]]>maxE) {
        maxE=energyBC[BCs[ibc]];
        seed=BCs[ibc];
      }
    }
    return seed;
  } else {
    for(int ibc=0;ibc<nPFBC;++ibc) {
      if(indexSCPFBC[ibc]==sc) BCs.push_back(ibc);
    }
    // find the seed BC
    // N.B. Since the PFBC collection was set to be the one made by PFphoton reco, then sometimes there could be not the seed of the electron PFSC. to be fixed
    float maxE=-1.;
    int seed=-1;
    for(unsigned int ibc=0; ibc<BCs.size(); ++ibc) { 
      if(energyPFBC[BCs[ibc]]>maxE) {
        maxE=energyPFBC[BCs[ibc]];
        seed=BCs[ibc];
      }
    }
    return seed;
  }
}
