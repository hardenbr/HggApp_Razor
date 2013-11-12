// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>

#include "include/HZZEleIDSelector.hh"
#include "include/eIDCiChzzSelector.hh"
#include "include/eIDSimpleCutsSelector.hh"
#include "include/EGammaMvaEleEstimator.h"
#include "include/ElectronEffectiveArea.h"
#include "../HiggsAnalysisTools/macro/LumiReweightingStandAlone.h"

using namespace std;
using namespace reweight;

enum idType {
  kIsoHWW2011 = 0, // HWW cuts 2011
  kBDTHWW2011_withIP, // HWW cuts 2011
  kIso,
  kIsoEACorr,
  kBDTHZZ_withIP
};

enum cutBasedIdType {
  kVeto = 0,
  kLoose,
  kMedium,
  kTight
};
  

float Aeff_neu_dr04[7], Aeff_pho_dr04[7];

float Aeff_neu_dr04_2011[7] = { 0.045, 0.065, 0.068, 0.057, 0.058, 0.061, 0.110 };
float Aeff_pho_dr04_2011[7] = { 0.140, 0.130, 0.079, 0.130, 0.150, 0.160, 0.180 };

float Aeff_neu_dr04_2012[7] = { 0.041, 0.068, 0.075, 0.068, 0.071, 0.078, 0.140 };
float Aeff_pho_dr04_2012[7] = { 0.144, 0.138, 0.084, 0.155, 0.201, 0.223, 0.265 };

// H->ZZ detector based effective areas
float Aeff_ecal_dr03[2] = { 0.101, 0.046 };
float Aeff_hcal_dr03[2] = { 0.021, 0.040 };

bool cicidval(int *cic, int level) { return (cic[level]>>0)%2; }
bool cicisoval(int *cic, int level) { return (cic[level]>>1)%2; }
bool cicconvval(int *cic, int level) { return (cic[level]>>2)%2; }
bool cicipval(int *cic, int level) { return (cic[level]>>3)%2; }

bool passHWWID(Float_t eta, Float_t pt, Float_t bdthww, Float_t bdthzz, Float_t rho, Float_t combIso, Float_t combPFIsoHWW, idType type) {
  if(type == kIsoHWW2011) {
    if(fabs(eta)<1.479) return (combPFIsoHWW/pt < 0.13);
    else return (combPFIsoHWW/pt < 0.09);
  }

  if(type == kIsoEACorr) {
    if(pt>20) {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.23); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.20);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.12);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.11);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return (combIso/pt < 0.049);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < 0.070);
      if(fabs(eta) >=  2.4) return (combIso/pt < 0.010);
    } else {
      if(fabs(eta) <  1.0) return (combIso/pt < 0.20); 
      if(fabs(eta) >=  1.0 && fabs(eta) < 1.479) return (combIso/pt < 0.21);
      if(fabs(eta) >=  1.479 && fabs(eta) < 2.0) return (combIso/pt < 0.13);
      if(fabs(eta) >=  2.0 && fabs(eta) < 2.2) return (combIso/pt < 0.10);
      if(fabs(eta) >=  2.2 && fabs(eta) < 2.3) return(combIso/pt < -0.04);
      if(fabs(eta) >=  2.3 && fabs(eta) < 2.4) return (combIso/pt < -0.03);
      if(fabs(eta) >=  2.4) return (combIso/pt < -0.03);
    }
  }

  if(type == kIso) {
    if(fabs(eta) < 1.479) return (combIso/pt < 0.29);
    else return (combIso/pt < 0.21);
  }

  if(type == kBDTHWW2011_withIP) {
    if(pt < 20 && fabs(eta) < 1.0) return (bdthww > 0.139);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.525);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.543);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthww > 0.947);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthww > 0.950);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthww > 0.884);
  }

  if(type == kBDTHZZ_withIP) {
    // WP with same fake rate as HWW with IP
    if(pt < 20 && fabs(eta) < 1.0) return (bdthzz > 0.075);
    if(pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.075);
    if(pt < 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.091);
    if(pt >= 20 && fabs(eta) < 1.0) return (bdthzz > 0.064);
    if(pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) return (bdthzz > 0.071);
    if(pt >= 20 && fabs(eta) >= 1.479 && fabs(eta) < 2.500) return (bdthzz > 0.067);
  }

  return false;
}

float isoEAWW2012(float eta, float pfIsoChHad04, float pfIsoNHad04_NoEA, float pfIsoPhoton04_NoEA, float rho) {
  float abseta=fabs(eta);

  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = ElectronEffectiveArea::kEleEAData2012;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGammaAndNeutralHad_ =  ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04;

  float eff_area_ga  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
  float eff_area_nh  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
  float eff_area_ganh = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGammaAndNeutralHad_, abseta, effAreaTarget_);

  float iso = pfIsoChHad04;
  iso += max<float>(0.,pfIsoNHad04_NoEA+pfIsoPhoton04_NoEA - eff_area_ganh*rho);
  return iso;
}

bool passHWWEleId2012(float pt, float eta, float bdt, float isoEA, bool matchConv, float misshits, float d0, float dz, int step) {
  // steps: 0=id, 1=iso, 2=conv, 3=ip, 4=full
  float abseta=fabs(eta);
  bool ELE_ID_EGAMMA_2012 = (
			     ( pt <= 20 && abseta >= 0.000 && abseta < 0.800 && bdt > 0.00 ) ||
                             ( pt <= 20 && abseta >= 0.800 && abseta < 1.479 && bdt > 0.10 ) ||
                             ( pt <= 20 && abseta >= 1.479 && abseta < 2.500 && bdt > 0.62 ) ||
                             ( pt >  20 && abseta >= 0.000 && abseta < 0.800 && bdt > 0.94 ) ||
                             ( pt >  20 && abseta >= 0.800 && abseta < 1.479 && bdt > 0.85 ) ||
                             ( pt >  20 && abseta >= 1.479 && abseta < 2.500 && bdt > 0.92 )
                             );
  
  bool ELE_ISO_EGAMMA_2012 = (isoEA/pt < 0.15);
  bool ELE_CONV_2012 = (!matchConv && misshits==0);
  bool ELE_IP_2012 = (fabs(d0)<0.02 && fabs(dz)<0.1);
  
  if(step==0) return ELE_ID_EGAMMA_2012;
  if(step==1) return ELE_ISO_EGAMMA_2012;
  if(step==2) return ELE_CONV_2012;
  if(step==3) return ELE_IP_2012;
  if(step==4) return (ELE_ID_EGAMMA_2012 && ELE_ISO_EGAMMA_2012 && ELE_CONV_2012 && ELE_IP_2012);
}

bool passHZZ4lEleId2012(float pt, float eta, float bdt, float isoEA, float misshits, float sip, int step) {
  // steps: 0=id, 1=iso, 2=conv, 3=ip, 4=full
  float abseta=fabs(eta);
  bool ELE_ID_EGAMMA_2012 = (
			     ( pt <= 10 && abseta >= 0.000 && abseta < 0.800 && bdt > 0.470 ) ||
                             ( pt <= 10 && abseta >= 0.800 && abseta < 1.479 && bdt > 0.004 ) ||
                             ( pt <= 10 && abseta >= 1.479 && abseta < 2.500 && bdt > 0.295 ) ||
                             ( pt >  10 && abseta >= 0.000 && abseta < 0.800 && bdt > 0.500 ) ||
                             ( pt >  10 && abseta >= 0.800 && abseta < 1.479 && bdt > 0.120 ) ||
                             ( pt >  10 && abseta >= 1.479 && abseta < 2.500 && bdt > 0.600 )
                             );
  
  bool ELE_ISO_EGAMMA_2012 = (isoEA/pt < 0.40);
  bool ELE_CONV_2012 = (misshits<=1);
  bool ELE_IP_2012 = (fabs(sip)<4);
  
  if(step==0) return ELE_ID_EGAMMA_2012;
  if(step==1) return ELE_ISO_EGAMMA_2012;
  if(step==2) return ELE_CONV_2012;
  if(step==3) return ELE_IP_2012;
  if(step==4) return (ELE_ID_EGAMMA_2012 && ELE_ISO_EGAMMA_2012 && ELE_CONV_2012 && ELE_IP_2012);
}

void makeFriendHZZIsolation(const char* file, int ismc) {

  // 2012
  LumiReWeighting LumiWeights( "data/s10MCPileUp.root",
                               "data/PUdata2012Final.root",
                               "pileup","pileup");

  // 2011
  //   LumiReWeighting LumiWeights( "/afs/cern.ch/user/e/emanuele/workspace/public/pileup/s6MCPileUp.root",
  //                                "/afs/cern.ch/user/e/emanuele/workspace/public/pileup/PUTarget.Full2011.160404-180252.root",
  //                                "pileup","pileup");

  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = ElectronEffectiveArea::kEleEAData2011;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGammaAndNeutralHad_ =  ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04;

  TFile *pF = TFile::Open(file);
  TTree *pT = (TTree*)pF->Get("eleIDdir/T1");
  
  TString nF(file);
  nF.ReplaceAll(".root","_hzzisoFriend.root");
  TFile *fF = TFile::Open(nF,"recreate");

  EGammaMvaEleEstimator *fElectronIsoMVA = new EGammaMvaEleEstimator();
  vector<string> eleiso_weightfiles;
  eleiso_weightfiles.push_back("elebdtweights/ElectronIso_BDTG_V0_BarrelPt5To10.weights.xml");
  eleiso_weightfiles.push_back("elebdtweights/ElectronIso_BDTG_V0_EndcapPt5To10.weights.xml");
  eleiso_weightfiles.push_back("elebdtweights/ElectronIso_BDTG_V0_BarrelPt10ToInf.weights.xml");
  eleiso_weightfiles.push_back("elebdtweights/ElectronIso_BDTG_V0_EndcapPt10ToInf.weights.xml");

  fElectronIsoMVA->initialize("EleIso_BDTG_IsoRings",
			      EGammaMvaEleEstimator::kIsoRings,
			      true,
			      eleiso_weightfiles);

  Float_t chaPFIso[8], neuPFIso[8], phoPFIso[8], rho, eta, pt;
  Float_t trkIso,ecalIso,hcalIso;
  Float_t npu[3];
  pT->SetBranchAddress("chaPFIso", chaPFIso);
  pT->SetBranchAddress("neuPFIso", neuPFIso);
  pT->SetBranchAddress("phoPFIso", phoPFIso);
  pT->SetBranchAddress("trkIso04", &trkIso);
  pT->SetBranchAddress("ecalIso04", &ecalIso);
  pT->SetBranchAddress("hcalIso04", &hcalIso);
  pT->SetBranchAddress("rho", &rho);
  pT->SetBranchAddress("eta", &eta);
  pT->SetBranchAddress("pt", &pt);
  pT->SetBranchAddress("npu", npu);

  fF->mkdir("eleIDdir");
  TTree *fT = new TTree("T1","tree with hzz isolation");
  Float_t combDetIso, combIso, combIsoNoEA;
  Float_t combIsoHww;
  Float_t mvaPFIso;
  Float_t puW;
  fT->Branch("combDetIsoHZZ",&combDetIso,"combDetIsoHZZ/F");
  fT->Branch("combPFIsoHZZ",&combIso,"combPFIsoHZZ/F");
  fT->Branch("combPFIsoHZZNoEA",&combIsoNoEA,"combPFIsoHZZNoEA/F");
  fT->Branch("mvaPFIso",&mvaPFIso,"mvaPFIso/F");
  fT->Branch("combIsoHww",&combIsoHww,"combIsoHww/F");
  fT->Branch("puW", &puW, "puW/F");

  for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);

     combIsoNoEA = chaPFIso[3]+neuPFIso[3]+phoPFIso[3]; // [0] = cone 0.1, then steps of 0.1 up to 0.7 [7] is directional isolation with DR=0.4

     float abseta=fabs(eta);
     float eff_area_ga  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
     float eff_area_nh  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
     float eff_area_ganh = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGammaAndNeutralHad_, abseta, effAreaTarget_);
     
     combIso = chaPFIso[3];
     combIso += max<float>(0.,neuPFIso[3]+phoPFIso[3] - eff_area_ganh*rho);

     int ieta = (fabs(eta)<1.479) ? 0 : 1;
     combDetIso = trkIso + ecalIso - rho * Aeff_ecal_dr03[ieta] + hcalIso - rho * Aeff_hcal_dr03[ieta];

     // slow and useless
     /*
     mvaPFIso = fElectronIsoMVA->isoMvaValue(pt,eta,rho,ElectronEffectiveArea::kEleEAData2011,
					     chaPFIso[0],chaPFIso[1]-chaPFIso[0],chaPFIso[2]-chaPFIso[1],chaPFIso[3]-chaPFIso[2],chaPFIso[4]-chaPFIso[3],
					     phoPFIso[0],phoPFIso[1]-phoPFIso[0],phoPFIso[2]-phoPFIso[1],phoPFIso[3]-phoPFIso[2],phoPFIso[4]-phoPFIso[3],
					     neuPFIso[0],neuPFIso[1]-neuPFIso[0],neuPFIso[2]-neuPFIso[1],neuPFIso[3]-neuPFIso[2],neuPFIso[4]-neuPFIso[3],
					     true);
     */
     mvaPFIso = -1.;

     combIsoHww = isoEAWW2012(eta,chaPFIso[3],neuPFIso[3],phoPFIso[3],rho);

     // this is only needed in MC
     if(ismc) puW = LumiWeights.weight(npu[1]);
     else puW = 1.0;

     fT->Fill();
  }

  fF->cd("eleIDdir");
  fT->Write();
  fF->Close();

  cout << "DONE. Friend tree is in file: " << nF.Data() << endl;

}

void makeFriendHZZIdBits(const char* file, int ismc) {


  TFile *pF = TFile::Open(file);
  TTree *pT = (TTree*)pF->Get("eleIDdir/T1");

  TString isofriend = TString(file);
  isofriend.ReplaceAll(".root","_hzzisoFriend.root");
  pT->AddFriend( "eleIDdir/isoT1 = eleIDdir/T1", isofriend );
 
  TString nF(file);
  nF.ReplaceAll(".root","_hzzidbitsFriend.root");
  TFile *fF = TFile::Open(nF,"recreate");

  Float_t eta, abseta, pt, rho, vertices;
  Float_t bdthww[2], newbdthww[4], combPFIsoHWW;
  Float_t combPFIsoHZZ, bdthzz[4];
  Float_t combIsoHww;
  Float_t chaPFIso[8], neuPFIso[8], phoPFIso[8], mvaPFIso;
  Float_t mass; // not dummy only for TP trees
  Int_t DenomFakeSmurf, ecalseed;
  Float_t eop,eseedopin,HoE,deta,dphi,see,fbrem,dist,dcot,d0,dz,sip,trkIso,ecalIso,hcalIso,IoEmIoP;
  Int_t missHits;
  Bool_t matchConv;
  Float_t npu[3];
  Int_t mcmatch;
  Float_t puW;
  pT->SetBranchAddress("bdthww", bdthww);
  pT->SetBranchAddress("newbdthww", newbdthww);
  pT->SetBranchAddress("bdthzz",bdthzz);
  pT->SetBranchAddress("combPFIsoHWW", &combPFIsoHWW);
  pT->SetBranchAddress("combPFIsoHZZ",&combPFIsoHZZ);
  pT->SetBranchAddress("chaPFIso", chaPFIso);
  pT->SetBranchAddress("neuPFIso", neuPFIso);
  pT->SetBranchAddress("phoPFIso", phoPFIso);
  pT->SetBranchAddress("DenomFakeSmurf", &DenomFakeSmurf);
  pT->SetBranchAddress("eta", &eta);
  pT->SetBranchAddress("pt", &pt);
  pT->SetBranchAddress("rho", &rho);
  pT->SetBranchAddress("vertices", &vertices);
  pT->SetBranchAddress("EoP", &eop);
  pT->SetBranchAddress("EseedoPin", &eseedopin);
  pT->SetBranchAddress("HoE", &HoE);
  pT->SetBranchAddress("deta", &deta);
  pT->SetBranchAddress("dphi", &dphi);
  pT->SetBranchAddress("see", &see);
  pT->SetBranchAddress("fbrem", &fbrem);
  pT->SetBranchAddress("IoEmIoP", &IoEmIoP);
  pT->SetBranchAddress("dist", &dist);
  pT->SetBranchAddress("dcot", &dcot);
  pT->SetBranchAddress("d0",&d0);
  pT->SetBranchAddress("dz",&dz);
  pT->SetBranchAddress("ip3ds",&sip);
  pT->SetBranchAddress("matchConv",&matchConv);
  pT->SetBranchAddress("missHits",&missHits);
  pT->SetBranchAddress("ecaldriven", &ecalseed);
  pT->SetBranchAddress("trkIso", &trkIso);
  pT->SetBranchAddress("ecalIso04", &ecalIso);
  pT->SetBranchAddress("hcalIso04", &hcalIso);
  pT->SetBranchAddress("mvaPFIso04", &mvaPFIso);
  pT->SetBranchAddress("combIsoHww", &combIsoHww);
  pT->SetBranchAddress("npu", npu);
  pT->SetBranchAddress("mcmatch", &mcmatch);
  pT->SetBranchAddress("puW", &puW);
  if(!TString(file).Contains("fake")) pT->SetBranchAddress("mass", &mass);
  else mass=-1.0;

  fF->mkdir("eleIDdir");
  TTree *fT = new TTree("T1","tree with hzz isolation");
  // the new WPs with full isolation for triggering electrons
  Int_t WP95trg, WP90trg, WP85trg, WP80trg, WP70trg;
  // the new WPs with charged only isolation for non triggering electrons
  Int_t WP95notrg, WP90notrg, WP85notrg, WP80notrg, WP70notrg;
  // the hww2011 WP and the one with the same efficiency
  // hzz WP is using the MVA for the unbiased electrons
  Int_t hwwWP, newhwwWP, newhzzWP;
  Int_t hwwWPisoonly, newhwwWPisoonly, newhzzWPisoonly;
  Int_t newhzzWPconvonly;
  Int_t hwwWPidonly, newhwwWPidonly, newhzzWPidonly;
  Int_t cicall[5], cicid[5], ciciso[5];
  Int_t spid[4];
  Int_t hzzMvaLoose, hzzMvaTight;
  // first 4 variables needed for TP
  fT->Branch("mass", &mass, "mass/F");
  fT->Branch("pt", &pt, "pt/F");
  fT->Branch("abseta", &abseta, "abseta/F");
  fT->Branch("vertices", &vertices, "vertices/F");
  fT->Branch("mcmatch", &mcmatch, "mcmatch/I");
  // for triggering eles
  fT->Branch("wp95trg", &WP95trg, "wp95trg/I");
  fT->Branch("wp90trg", &WP90trg, "wp90trg/I");
  fT->Branch("wp85trg", &WP85trg, "wp85trg/I");
  fT->Branch("wp80trg", &WP80trg, "wp80trg/I");
  fT->Branch("wp70trg", &WP70trg, "wp70trg/I");
  fT->Branch("newhwwWP", &newhwwWP, "newhwwWP/I"); // 2012 WP
  fT->Branch("newhwwWPisoonly", &newhwwWPisoonly, "newhwwWPisoonly/I"); // 2012 WP
  fT->Branch("newhwwWPidonly", &newhwwWPidonly, "newhwwWPidonly/I"); // 2012 WP
  fT->Branch("hwwWP", &hwwWP, "hwwWP/I");  // 2011 WP
  fT->Branch("hwwWPisoonly", &hwwWPisoonly, "hwwWPisoonly/I");  // 2011 WP
  fT->Branch("hwwWPidonly", &hwwWPidonly, "hwwWPidonly/I");  // 2011 WP
  // for non triggering eles
  fT->Branch("wp95notrg", &WP95notrg, "chwp95notrg/I");
  fT->Branch("wp90notrg", &WP90notrg, "chwp90notrg/I");
  fT->Branch("wp85notrg", &WP85notrg, "chwp85notrg/I");
  fT->Branch("wp80notrg", &WP80notrg, "chwp80notrg/I");
  fT->Branch("wp70notrg", &WP70notrg, "chwp70notrg/I");
  // same as HWW DenomFakeSmurf: change name for the friend tree
  fT->Branch("denom", &DenomFakeSmurf, "denom/I");
  // the cic used for HZZ
  fT->Branch("cicall", cicall, "cicall[5]/I");
  fT->Branch("cicid", cicid, "cicid[5]/I");
  fT->Branch("ciciso", ciciso, "ciciso[5]/I");
  fT->Branch("spid", spid, "spid[4]/I");
  fT->Branch("newhzzWP", &newhzzWP, "newhzzWP/I"); // 2012 WP
  fT->Branch("newhzzWPisoonly", &newhzzWPisoonly, "newhzzWPisoonly/I"); // 2012 WP
  fT->Branch("newhzzWPidonly", &newhzzWPidonly, "newhzzWPidonly/I"); // 2012 WP
  fT->Branch("newhzzWPconvonly", &newhzzWPconvonly, "newhzzWPconvonly/I"); // 2012 WP
  // mva 2012 duncan's WP
  fT->Branch("hzzMvaLoose", &hzzMvaLoose, "hzzMvaLoose/I");
  fT->Branch("hzzMvaTight", &hzzMvaTight, "hzzMvaTight/I");
  fT->Branch("puW", &puW, "puW/F");

  HZZEleIDSelector aSel;

  for(int i=0; i<pT->GetEntries(); i++) {
    if (i%10000 == 0) std::cout << ">>> Analyzing event # " << i << " / " << pT->GetEntries() << " entries" << std::endl;
     pT->GetEntry(i);

     abseta = fabs(eta);

     hwwWP=0;
     if(passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kBDTHWW2011_withIP) && 
	passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kIsoHWW2011)) hwwWP = 1;
     hwwWPisoonly=0;
     if(passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kIsoHWW2011)) hwwWPisoonly = 1;
     hwwWPidonly=0;
     if(passHWWID(eta,pt,bdthww[0],newbdthww[3],rho,combPFIsoHZZ,combPFIsoHWW,kBDTHWW2011_withIP)) hwwWPidonly = 1;

     WP95trg=WP90trg=WP85trg=WP80trg=WP70trg=0;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP95,HZZEleIDSelector::kMVABiased)) WP95trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP90,HZZEleIDSelector::kMVABiased)) WP90trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP85,HZZEleIDSelector::kMVABiased)) WP85trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVABiased)) WP80trg=1;
     if(aSel.output(pt,eta,newbdthww[3],combPFIsoHZZ,HZZEleIDSelector::kWP70,HZZEleIDSelector::kMVABiased)) WP70trg=1;

     // 2012 WP. Same FR as 2011
     newhwwWP=0;
     if(passHWWEleId2012(pt,eta,newbdthww[3],combIsoHww,matchConv,missHits,d0,dz,4)) newhwwWP=1;
     newhwwWPisoonly=0;
     if(passHWWEleId2012(pt,eta,newbdthww[3],combIsoHww,matchConv,missHits,d0,dz,1)) newhwwWPisoonly=1;
     newhwwWPidonly=0;
     if(passHWWEleId2012(pt,eta,newbdthww[3],combIsoHww,matchConv,missHits,d0,dz,0)) newhwwWPidonly=1;

     WP95notrg=WP90notrg=WP85notrg=WP80notrg=WP70notrg=0;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP95,HZZEleIDSelector::kMVAUnbiased)) WP95notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP90,HZZEleIDSelector::kMVAUnbiased)) WP90notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP85,HZZEleIDSelector::kMVAUnbiased)) WP85notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP80,HZZEleIDSelector::kMVAUnbiased)) WP80notrg=1;
     if(aSel.output(pt,eta,bdthzz[3],combPFIsoHZZ,HZZEleIDSelector::kWP70,HZZEleIDSelector::kMVAUnbiased)) WP70notrg=1;

     // CiC...
     eIDCiChzzSelector cicsel;
     int cic[5];
     int ieta = (fabs(eta)<1.479) ? 0 : 1;
     for(int i=0; i<5; i++) cic[i] = cicsel.ElectronId_V03(pt,eta,see,eop,eseedopin,fbrem,
							   trkIso,
							   ecalIso - rho * Aeff_ecal_dr03[ieta],
							   hcalIso - rho * Aeff_hcal_dr03[ieta],
							   d0,missHits,deta,dphi,HoE,dcot,
							   dist,!ecalseed,i,false);
     for(int i=0; i<5; i++) {
       cicall[i] = (cicidval(cic,i) && cicisoval(cic,i) && missHits<=1) ? 1 : 0;
       cicid[i] = (cicidval(cic,i) && missHits<=1) ? 1 : 0;
       ciciso[i] = (cicisoval(cic,i)) ? 1 : 0;
     }
     
     // simple cuts 2012
     eIDSimpleCutsSelector spsel(deta,dphi,see,HoE,d0,dz,IoEmIoP);
     for(int i=0; i<4; i++) {
       spid[i] = spsel.output(eta,i);
     }

     // H->ZZ 2012 cut
     newhzzWP=0;
     if(passHZZ4lEleId2012(pt,eta,bdthzz[3],combPFIsoHZZ,missHits,sip,4)) newhzzWP=1;
     newhzzWPisoonly=0;
     if(passHZZ4lEleId2012(pt,eta,bdthzz[3],combPFIsoHZZ,missHits,sip,1)) newhzzWPisoonly=1;
     newhzzWPidonly=0;
     if(passHZZ4lEleId2012(pt,eta,bdthzz[3],combPFIsoHZZ,missHits,sip,0)) newhzzWPidonly=1;
     newhzzWPconvonly=0;
     if(passHZZ4lEleId2012(pt,eta,bdthzz[3],combPFIsoHZZ,missHits,sip,2)) newhzzWPconvonly=1;

     // MVA iso rings and MVA ID
     hzzMvaLoose=hzzMvaTight=0;
     if(aSel.output(pt,eta,bdthzz[3],mvaPFIso,HZZEleIDSelector::kMVALoose,HZZEleIDSelector::kMVAUnbiased)) hzzMvaLoose=1;
     if(aSel.output(pt,eta,bdthzz[3],mvaPFIso,HZZEleIDSelector::kMVATight,HZZEleIDSelector::kMVAUnbiased)) hzzMvaTight=1;

     if(ismc==0) mcmatch=1;

     fT->Fill();
  }
  fF->cd("eleIDdir");
  fT->Write();
  fF->Close();

  cout << "DONE. Friend tree is in file: " << nF << endl;
}


int main(int argc, char* argv[]) {

  int year = 2011;
  for(int i=0; i<7; i++) {
    if(year==2011) {
      Aeff_neu_dr04[i] = Aeff_neu_dr04_2011[i];
      Aeff_pho_dr04[i] = Aeff_pho_dr04_2011[i];
    } else if(year==2012) {
      Aeff_neu_dr04[i] = Aeff_neu_dr04_2012[i];
      Aeff_pho_dr04[i] = Aeff_pho_dr04_2012[i];
    } else {
      cout << "wrong year set. Returning 0." << endl;
      return 0;
    }
  }


  char files1[500], files2[500], fileb1[500], fileb2[500], fileb3[500], fileb4[500], fileb5[500]; 
  sprintf(files1,"macro/results_data_2012/electrons.root");
  sprintf(files2,"macro/results_data_2012/electrons_zeemc.root");
  sprintf(fileb1,"macro/results_data_2012/fakes.root");
  sprintf(fileb2,"macro/results_data_2012/fakes-wlnu1e.root");
  sprintf(fileb3,"macro/results_data_2012/fakes-zll1e.root");
  sprintf(fileb4,"macro/results_data_2012/fakes-ewksub-wlnu.root");
  sprintf(fileb5,"macro/results_data_2012/fakes-ewksub-zee.root");

  cout << "\t===> DOING ISOLATION FRIEND TREES <===" << endl;
  // isolation
  makeFriendHZZIsolation(files1,0);
  makeFriendHZZIsolation(files2,1);
  makeFriendHZZIsolation(fileb1,0);
  makeFriendHZZIsolation(fileb2,0);
  makeFriendHZZIsolation(fileb3,0);
  makeFriendHZZIsolation(fileb4,1);
  makeFriendHZZIsolation(fileb5,1);

  cout << "\t===> DOING ID FRIEND TREES <===" << endl;
  // id bits
  makeFriendHZZIdBits(files1,0);
  makeFriendHZZIdBits(files2,1);
  makeFriendHZZIdBits(fileb1,0);
  makeFriendHZZIdBits(fileb2,0);
  makeFriendHZZIdBits(fileb3,0);
  makeFriendHZZIdBits(fileb4,1);
  makeFriendHZZIdBits(fileb5,1);
}
