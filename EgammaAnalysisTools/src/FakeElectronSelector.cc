
#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/FakeElectronSelector.hh"
#include "CommonTools/include/Counters.hh"
#include "EgammaAnalysisTools/include/eIDCiChzzSelector.hh"

using namespace bits;
using namespace std;

FakeElectronSelector::FakeElectronSelector(TTree *tree)
  : Egamma(tree) {
  
  _isData = true;       // chiara
  
  // configuring the electron BDT
  fMVAHWW = new ElectronIDMVA();
  fMVAHWW->Initialize("BDTG method",
                      "elebdtweights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml" ,                
                      ElectronIDMVA::kWithIPInfo);

  fMVAHWWWithIso = new ElectronIDMVA();
  fMVAHWWWithIso->Initialize("BDTG method",
			     "elebdtweights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml" ,                
			     ElectronIDMVA::kIDIsoCombined);
  
  // configuring the electron BDT for H->ZZ
  fMVAHZZDanV0 = new ElectronIDMVAHZZ();
  fMVAHZZSiV0 = new ElectronIDMVAHZZ();
  fMVAHZZSiV1 = new ElectronIDMVAHZZ();
  fMVAHZZSiDanV2 = new ElectronIDMVAHZZ();

  // New H->ZZ unbiased DATA training, Daniele's variables
  fMVAHZZDanV0->Initialize("BDTSimpleCat",
                           "elebdtweights/DanieleMVA_BDTCat_BDTG_DanV0.weights.xml",
                           ElectronIDMVAHZZ::kBDTDanV0);

  // New H->ZZ unbiased DATA training, Si's HWW 2011 variables
  fMVAHZZSiV0->Initialize("BDTSimpleCat",
                          "elebdtweights/DanieleMVA_BDTCat_BDTG_SiV0.weights.xml",
                          ElectronIDMVAHZZ::kBDTSiV0);

  // New H->ZZ unbiased DATA training, Si's HWW 2012 variables
  fMVAHZZSiV1->Initialize("BDTSimpleCat",
                          "elebdtweights/DanieleMVA_BDTCat_BDTG_SiV1.weights.xml",
                          ElectronIDMVAHZZ::kBDTSiV1);

  // New H->ZZ unbiased DATA training, Daniele's + Si's variables 
  fMVAHZZSiDanV2->Initialize("BDTSimpleCat",
                             "elebdtweights/DanieleMVA_BDTCat_BDTG_SiDanV2.weights.xml",
                             ElectronIDMVAHZZ::kBDTSiDanV2);


  // configuring the electron BDT for H->WW
  fMVAHWWDanV0 = new ElectronIDMVAHZZ();
  fMVAHWWSiV0 = new ElectronIDMVAHZZ();
  fMVAHWWSiV1 = new ElectronIDMVAHZZ();
  fMVAHWWSiDanV2 = new ElectronIDMVAHZZ();

  // New H->WW unbiased DATA training, Daniele's variables
  fMVAHWWDanV0->Initialize("BDTSimpleCat",
                           "elebdtweights/DanieleMVA_DenomHWW_BDTCat_BDTG_DanV0.weights.xml",
                           ElectronIDMVAHZZ::kBDTHWWDanV0);

  // New H->WW unbiased DATA training, Si's HWW 2011 variables
  fMVAHWWSiV0->Initialize("BDTSimpleCat",
                          "elebdtweights/DanieleMVA_DenomHWW_BDTCat_BDTG_SiV0.weights.xml",
                          ElectronIDMVAHZZ::kBDTHWWSiV0);

  // New H->WW unbiased DATA training, Si's HWW 2012 variables
  fMVAHWWSiV1->Initialize("BDTSimpleCat",
                          "elebdtweights/DanieleMVA_DenomHWW_BDTCat_BDTG_SiV1.weights.xml",
                          ElectronIDMVAHZZ::kBDTHWWSiV1);

  // New H->WW unbiased DATA training, Daniele's + Si's variables 
  fMVAHWWSiDanV2->Initialize("BDTSimpleCat",
                             "elebdtweights/DanieleMVA_DenomHWW_BDTCat_BDTG_SiDanV2.weights.xml",
                             ElectronIDMVAHZZ::kBDTHWWSiDanV2);


  // chiara
  // to read good run list
  if (_isData) {
    std::string goodRunGiasoneFile = "config/json/goodCollisions2012.json";
    setJsonGoodRunList(goodRunGiasoneFile); 
    fillRunLSMap();
  }
  
  // counter initialize
  myCounter.SetTitle("EVENT_COUNTER");
  myCounter.AddVar("event");
  myCounter.AddVar("trigger");
  myCounter.AddVar("denom");
  myCounter.AddVar("met");
  myCounter.AddVar("trasvMass");
  myCounter.AddVar("Zmass");
  myCounter.AddVar("leadingExist");
  myCounter.AddVar("deltaPhi");
  myCounter.AddVar("leadingPT");
}

FakeElectronSelector::~FakeElectronSelector() { }


void FakeElectronSelector::Loop(const char *outname) {

  // study vs eta
  float minEta = -2.5;
  float maxEta =  2.5;

  // ---------------------------------------------------------------------
  // to study the event selection
  char filename[200];
  sprintf(filename,"%s_FakeKineTree.root",outname);
  myOutKineTree = new FakeTree(filename);
  sprintf(filename,"%s_FakeIDTree.root",outname);
  myOutIDTree = new RedEleIDTree(filename);
  myOutIDTree->addElectronIdBits();
  myOutIDTree->addDenominatorFakeBits();
  myOutIDTree->addIsolations();
  myOutIDTree->addRunInfos();
  myOutIDTree->addMore();
  myOutIDTree->addTrackMomenta();

  // trigger: electrons - chiara
  cout << "using electrons triggers" << endl;
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v1");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v2");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v3");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v4");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v5");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v6");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v7");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v8");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v9");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v10");
  // 
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v1");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v2");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v3");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v4");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v5");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v6");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v7");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v8");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v9");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v10");
  // 
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v5");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v6");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v7");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v8");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v9");
  // 
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v1");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v2");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v3");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v4");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v5");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v6");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v7");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v8");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v9");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v10");
  // 
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v1");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v2");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v3");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v4");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v5");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v6");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v7");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v8");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v9");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIDL_TrkIsoVL_v10");

  // 2012 triggers
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_TrkIdVL_v");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  requiredTriggers.push_back("HLT_Ele8_CaloIdT_TrkIdVL_Jet30_v");
  requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  requiredTriggers.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v");

  // loop on events
  unsigned int lastLumi = 0;
  unsigned int lastRun  = 0;
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // good runs selection 
    if (_isData && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun  = runNumber;
        lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (_isData && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    // skipping even numbered events (2011A, used for eleID BDT training) - chiara
    // if (runNumber<=173692 && eventNumber%2==0) continue;
    
    // all events
    myCounter.IncrVar("event",1);

    // event selection: trigger
    reloadTriggerMask(true);           
    bool passedHLT = hasPassedHLT();   
    if ( _isData && !passedHLT ) continue;
    myCounter.IncrVar("trigger",1);

    // electrons passing the denominator selection (removed to get it unbiased. Set a bit at the end)
    std::vector<int> denomElectrons;
    for(int iele=0; iele<nEle; iele++) {
      // bool isGoodDenom = isDenomFake_smurfs(iele);             
      bool isGoodDenom = true; // passthrough and set the bit
      if (!isGoodDenom) continue;
      denomElectrons.push_back(iele);
    }
    
    // event selection: at least one candidate for denominator
    if (denomElectrons.size()==0) continue;
    myCounter.IncrVar("denom",1);    

    // denominators: kinematics
    std::pair<int,int> possibleDenom = getBestGoodElePair(denomElectrons);     
    int theDenom1(possibleDenom.first);
    int theDenom2(possibleDenom.second);
    TLorentzVector tlvDenom1, tlvDenom2;
    tlvDenom1.SetXYZT(pxEle[theDenom1],pyEle[theDenom1],pzEle[theDenom1],energyEle[theDenom1]);
    tlvDenom2.SetXYZT(pxEle[theDenom2],pyEle[theDenom2],pzEle[theDenom2],energyEle[theDenom2]);
    TVector3 tv3Denom1, tv3Denom2;
    tv3Denom1.SetXYZ(pxEle[theDenom1],pyEle[theDenom1],pzEle[theDenom1]);
    tv3Denom2.SetXYZ(pxEle[theDenom2],pyEle[theDenom2],pzEle[theDenom2]);
    if (theDenom1==-1) cout << "sanity check: impossibile!" << endl;

    // look for the leading jet not matching the HLT object - to further reduce the W contribution
    float maxEt = -1.;
    int leadingJet = -1;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      // bool HLTmatch = triggerMatch(p3Jet.Eta(),p3Jet.Phi(),0.2);      
      // if (HLTmatch) continue;                                         
      if ( fabs(p3Jet.Eta()) < maxEta && p3Jet.Pt() > maxEt) {   
	maxEt = p3Jet.Pt();
	leadingJet = jet;
      }
    }

    // variables used in the selection
    TVector3 p3Met(pxPFMet[0],pyPFMet[0],0.0);
    TVector3 p3LeadingCut(pxAK5PFPUcorrJet[leadingJet],pyAK5PFPUcorrJet[leadingJet],pzAK5PFPUcorrJet[leadingJet]);
    float WmT = sqrt(2*tlvDenom1.Pt()*p3Met.Pt()*(1-cos(tv3Denom1.Angle(p3Met))) );      
    float theInvMass = -10.;
    if (theDenom1>-1 && theDenom2>-1) theInvMass = (tlvDenom1+tlvDenom2).M();
    float deltaPhi = fabs(p3LeadingCut.DeltaPhi(tv3Denom1));

    // fill the kine tree - after HLT and denominator
    myOutKineTree -> fill( p3Met.Pt(), WmT, theInvMass, maxEt, tv3Denom1.Pt(), deltaPhi);
    myOutKineTree -> store();
    
    // event selection: met cut to reduce W->enu
    if( p3Met.Pt() > 20 ) continue;       
    myCounter.IncrVar("met",1);    

    // event selection: mT cut to reduce W->enu [highest pT denominator candidate used] // chiara
    // if (WmT > 25. ) continue;  
    myCounter.IncrVar("trasvMass",1);    

    // event selection: invariant mass between two electrons passing the denominator selection to reduce Z->ee and low resonances
    if (theInvMass>60. && theInvMass<120.) continue;  
    if (theInvMass>0.1 && theInvMass<12.)  continue;  
    myCounter.IncrVar("Zmass",1);    

    /* remove the jet pt cuts 
    // need at least one reco jet
    if ( leadingJet < 0 ) continue;
    myCounter.IncrVar("leadingExist",1);    

    // the leading jet must be back-to-back wrt the fake candidate
    if ( deltaPhi<1 ) continue;
    myCounter.IncrVar("deltaPhi",1);  
    
    // minimal cut on leading jet or photon ET
    // if (p3LeadingCut.Pt()<15) continue;   
    // if (p3LeadingCut.Pt()<20) continue;   
    // if (p3LeadingCut.Pt()<30) continue;        // chiara
    if (p3LeadingCut.Pt()<35) continue;   
    // if (p3LeadingCut.Pt()<50) continue;   
    myCounter.IncrVar("leadingPT",1);    
    */
    float ptLeadingJet = -1.0;
    if(leadingJet >= 0) ptLeadingJet = p3LeadingCut.Pt();

    // consider as denominator all the reco electrons matching the HLT candidate 
    // passing some requirements on the following variables:
    // Gsf Electron (Ecal or tracker driven)
    // - ecal, hcal and tracker Isolation as in the trigger
    // - H/E, sigma ieta ieta, deltaEta, deltaPhi
    // - full IP
    // - full conversion rejection 

    // weight to the events to reweight for jet pT
    // float theWeight = 1.;
    
    // fill the denominator: take only the highest pT denominator candidate
    Utils anaUtils;
    float etaFake = fabs(tv3Denom1.Eta()); 
    float phiFake  = tv3Denom1.Phi();
    float etFake  = tv3Denom1.Pt();

    // for the systematics of HWW - init
    float dREleJet_min = 1000;
    int closestJet=-1;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.Eta()) < 2.5 && p3Jet.Pt() > 10.0 ) {          
	float dREleJet = p3Jet.DeltaR(tv3Denom1);
	if(dREleJet<dREleJet_min) {
	  closestJet=jet;
	  dREleJet_min=dREleJet;
	}}}
    // for the systematics of HWW - end

    // does this denominator pass the IP cut as for H->WW ?  // chiara, hardcoded
    bool isDenomIP = true;
    int gsfTrack = gsfTrackIndexEle[theDenom1];
    float dxyEle = transvImpactParGsfTrack[gsfTrack];
    float dzEle  = eleDzPV(theDenom1,0);
    if ( fabs(dxyEle)>0.02 ) isDenomIP = false;
    if ( fabs(dzEle)>0.10 )  isDenomIP = false;
    // does this denominator pass the IP cut as for H->WW ?  

    // some eleID variables
    float HoE, s1s9, s9s25, phiwidth, etawidth, deta, dphi, fbrem, see, spp, eleopout, eopout, eop, eseedopin, nbrems, recoFlag, EleSCEta, EleSCPhi;
    float EleSCEt,EtaSeed,PhiSeed,ESeed,IEtaSeed,IPhiSeed,EtaCrySeed,PhiCrySeed,IEtaCrySeed,IPhiCrySeed;
    float oneoveremoneoverp, eledeta, d0, ip3d, ip3ds, kfhits, kflayers, kfchi2, e1x5e5x5, dcot, dist;
    float detacalo, dphicalo, sep, dz, gsfchi2, emax, etop, ebottom, eleft, eright,
      e2nd, e2x5right, e2x5left, e2x5top, e2x5bottom, 
      e2x5max, e1x5, e2x2, e3x3, e5x5, r9, nclu,
      phi, scenergy, scrawenergy, scesenergy;

    int kfTrack = trackIndexEle[theDenom1];

    // different p estimations
    float pcomb=tv3Denom1.Mag();
    TVector3 p3ModeGsf(pxModeGsfTrack[gsfTrack],pyModeGsfTrack[gsfTrack],pzModeGsfTrack[gsfTrack]);
    float pmodegsf=p3ModeGsf.Mag();
    TVector3 p3MeanGsf(pxGsfTrack[gsfTrack],pyGsfTrack[gsfTrack],pzGsfTrack[gsfTrack]);
    float pmeangsf=p3MeanGsf.Mag();
    TVector3 p3MeanKf(pxTrack[kfTrack],pyTrack[kfTrack],pzTrack[kfTrack]);
    float pmeankf=p3MeanKf.Mag();
    float pterrorgsf = ptErrorGsfTrack[gsfTrack];
    float pterrorkf = (kfTrack>-1) ? ptErrorTrack[kfTrack] : -1.0;
    float perrorele = trackMomentumErrorEle[theDenom1];

    double gsfsign   = (-eleDxyPV(theDenom1,0) >=0 ) ? 1. : -1.;
    bool matchConv = hasMatchedConversionEle[theDenom1];

    d0 = gsfsign * transvImpactParGsfTrack[gsfTrack];
    dz = eleDzPV(theDenom1,0);
    ip3d = gsfsign * impactPar3DGsfTrack[gsfTrack];
    ip3ds = ip3d/impactPar3DErrorGsfTrack[gsfTrack];
    kfchi2 = (kfTrack>-1) ? trackNormalizedChi2Track[kfTrack] : 0.0;
    kflayers = (kfTrack>-1) ? trackerLayersWithMeasurementTrack[kfTrack] : -1.0;
    kfhits = (kfTrack>-1) ? trackValidHitsTrack[kfTrack] : -1.0;
    gsfchi2 = trackNormalizedChi2GsfTrack[gsfTrack];
    int misshits = expInnerLayersGsfTrack[gsfTrack];
    dcot = convDistEle[theDenom1];
    dist = convDcotEle[theDenom1];
    bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[theDenom1], isEcalDriven);
    int ecalseed;
    HoE = hOverEEle[theDenom1];
    eledeta = deltaEtaEleClusterTrackAtCaloEle[theDenom1];
    deta = deltaEtaAtVtxEle[theDenom1];
    dphi = deltaPhiAtVtxEle[theDenom1];
    detacalo = deltaEtaAtCaloEle[theDenom1];
    dphicalo = deltaPhiAtCaloEle[theDenom1];
    fbrem = fbremEle[theDenom1];
    nbrems = nbremsEle[theDenom1];
    eleopout = eEleClusterOverPoutEle[theDenom1];
    eopout = eSeedOverPoutEle[theDenom1];
    eop = eSuperClusterOverPEle[theDenom1];
    if(ecaldriven) {
      ecalseed = 1;
      int sc = superClusterIndexEle[theDenom1];
      float seedEnergy = seedClusterEnergySC[sc];
      nclu = float(nBCSC[sc]);
      s1s9 = eMaxSC[sc]/eMaxSC[sc];
      s9s25 = e3x3SC[sc]/e5x5SC[sc];
      e1x5e5x5 = (e5x5SC[sc] - e1x5SC[sc])/e5x5SC[sc];
      phiwidth = phiWidthSC[sc];
      etawidth = etaWidthSC[sc];
      see = sqrt(covIEtaIEtaSC[sc]);
      sep = covIEtaIPhiSC[sc]/(sqrt(covIEtaIEtaSC[sc])*sqrt(covIPhiIPhiSC[sc]));
      spp = sqrt(covIPhiIPhiSC[sc]);
      oneoveremoneoverp = 1./energySC[sc]  - 1./tv3Denom1.Mag();
      eseedopin = seedEnergy/p3ModeGsf.Mag();
      emax = eMaxSC[sc];
      etop = eTopSC[sc];
      ebottom = eBottomSC[sc];
      eleft = eLeftSC[sc];
      eright = eRightSC[sc];
      e2nd = e2ndSC[sc];
      e2x5right = e2x5RightSC[sc];
      e2x5left = e2x5LeftSC[sc];
      e2x5top = e2x5TopSC[sc];
      e2x5bottom = e2x5BottomSC[sc];
      e2x5max = e2x5MaxSC[sc];
      e1x5 = e1x5SC[sc];
      e2x2 = e2x2SC[sc];
      e3x3 = e3x3SC[sc];
      e5x5 = e5x5SC[sc];
      r9 = e3x3SC[sc]/rawEnergySC[sc];
      recoFlag = recoFlagSC[sc];
      EleSCEta = etaSC[sc];
      EleSCPhi = phiSC[sc];            
      EleSCEt = energySC[sc]*fabs(sin(thetaSC[sc]));
      scenergy = energySC[sc];
      scrawenergy = rawEnergySC[sc];
      scesenergy = esEnergySC[sc];
      int seedclu = indexSeedBC(sc,ecaldriven);
      EtaSeed=etaBC[seedclu];
      PhiSeed=phiBC[seedclu];
      ESeed=energyBC[seedclu];
      IEtaSeed=iEtaBC[seedclu];
      IPhiSeed=iPhiBC[seedclu];
      EtaCrySeed=etaCrystalBC[seedclu];
      PhiCrySeed=phiCrystalBC[seedclu];
      IEtaCrySeed=seedXSC[sc];
      IPhiCrySeed=seedYSC[sc];
    } else {
      ecalseed = 0;
      int sc = PFsuperClusterIndexEle[theDenom1];
      if(sc>-1) {
	float seedEnergy = seedEnergyPFSC[sc];
        nclu = float(nBCPFSC[sc]);
        s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
        s1s9 = eMaxPFSC[sc]/eMaxPFSC[sc];
        e1x5e5x5 = (e5x5PFSC[sc] - e1x5PFSC[sc])/e5x5PFSC[sc];
        phiwidth = phiWidthPFSC[sc];
        etawidth = etaWidthPFSC[sc];
        see = sqrt(covIEtaIEtaPFSC[sc]);
        sep = covIEtaIPhiPFSC[sc]/(sqrt(covIEtaIEtaPFSC[sc])*sqrt(covIPhiIPhiPFSC[sc]));
        spp = sqrt(covIPhiIPhiPFSC[sc]);
        oneoveremoneoverp = 1./energyPFSC[sc]  - 1./tv3Denom1.Mag();
	eseedopin = seedEnergy/p3ModeGsf.Mag();
        emax = eMaxPFSC[sc];
        etop = eTopPFSC[sc];
        ebottom = eBottomPFSC[sc];
        eleft = eLeftPFSC[sc];
        eright = eRightPFSC[sc];
        e2nd = e2ndPFSC[sc];
        e2x5right = e2x5RightPFSC[sc];
        e2x5left = e2x5LeftPFSC[sc];
        e2x5top = e2x5TopPFSC[sc];
        e2x5bottom = e2x5BottomPFSC[sc];
        e2x5max = e2x5MaxPFSC[sc];
        e1x5 = e1x5PFSC[sc];
        e2x2 = e2x2PFSC[sc];
        e3x3 = e3x3PFSC[sc];
        e5x5 = e5x5PFSC[sc];
        r9 = e3x3PFSC[sc]/rawEnergyPFSC[sc];
        recoFlag = recoFlagPFSC[sc];
        EleSCEta = etaPFSC[sc];
        EleSCPhi = phiPFSC[sc];     
        EleSCEt = energyPFSC[sc]*fabs(sin(thetaPFSC[sc]));
        scenergy = energyPFSC[sc];
        scrawenergy = rawEnergyPFSC[sc];
        scesenergy = esEnergyPFSC[sc];
        int seedclu = indexSeedBC(sc,ecaldriven);
        EtaSeed=etaPFBC[seedclu];
        PhiSeed=phiPFBC[seedclu];
        ESeed=energyPFBC[seedclu];
        IEtaSeed=iEtaPFBC[seedclu];
        IPhiSeed=iPhiPFBC[seedclu];
        EtaCrySeed=etaCrystalPFBC[seedclu];
        PhiCrySeed=phiCrystalPFBC[seedclu];
        IEtaCrySeed=seedXPFSC[sc];
        IPhiCrySeed=seedYPFSC[sc];
      } else {
        s9s25 = 999.;
        see = 999.;
        spp = 999.;
      }
    }

    // this only happens in 42X EE because of a bug in the EcalClusterTools
    if(isnan(see) || isnan(spp) || isnan(sep)) continue;

    float pt = tlvDenom1.Pt();
    float eta = tlvDenom1.Eta();
    int charge = chargeEle[theDenom1];

    // CiC...
    eIDCiChzzSelector cicsel;
    int cic[5];
    int iecal = (fabs(eta)<1.479) ? 0 : 1;
    for(int i=0; i<5; i++) cic[i] = cicsel.ElectronId_V03(pt,EleSCEta,see,eop,eseedopin,fbrem,
							  dr03TkSumPtEle[theDenom1],
							  dr03EcalRecHitSumEtEle[theDenom1] - rhoFastjet * Aeff_ecal_dr03[iecal],
							  dr03HcalTowerSumEtEle[theDenom1] - rhoFastjet * Aeff_hcal_dr03[iecal],
							  d0,misshits,deta,dphi,HoE,dcot,
							  dist,!ecalseed,i,false);

    // some MVAs...
    float pfmva = pflowMVAEle[theDenom1];
    float lh=eleIdLikelihoodEle[theDenom1];
    float hwwbdts[2];
    hwwbdts[0] = eleBDT(fMVAHWW,theDenom1);
    hwwbdts[1] = eleBDTWithIso(fMVAHWWWithIso,theDenom1);
    float hzzbdts[4];
    hzzbdts[0] = eleBDT(fMVAHZZDanV0,theDenom1);
    hzzbdts[1] = eleBDT(fMVAHZZSiV0,theDenom1);
    hzzbdts[2] = eleBDT(fMVAHZZSiV1,theDenom1);
    //hzzbdts[3] = eleBDT(fMVAHZZSiDanV2,theDenom1);
    hzzbdts[3] = mvaidnontrigEle[theDenom1];
    float newhwwbdts[4];
    newhwwbdts[0] = eleBDT(fMVAHWWDanV0,theDenom1);
    newhwwbdts[1] = eleBDT(fMVAHWWSiV0,theDenom1);
    newhwwbdts[2] = eleBDT(fMVAHWWSiV1,theDenom1);
    //newhwwbdts[3] = eleBDT(fMVAHWWSiDanV2,theDenom1);
    newhwwbdts[3] = mvaidtrigEle[theDenom1];

    // isolations
    float chaPfIso[8], phoPfIso[8], neuPfIso[8];
    chaPfIso[0]=pfCandChargedIso01Ele[theDenom1];
    phoPfIso[0]=pfCandPhotonIso01Ele[theDenom1];
    neuPfIso[0]=pfCandNeutralIso01Ele[theDenom1];
    
    chaPfIso[1]=pfCandChargedIso02Ele[theDenom1];
    phoPfIso[1]=pfCandPhotonIso02Ele[theDenom1];
    neuPfIso[1]=pfCandNeutralIso02Ele[theDenom1];
    
    chaPfIso[2]=pfCandChargedIso03Ele[theDenom1];
    phoPfIso[2]=pfCandPhotonIso03Ele[theDenom1];
    neuPfIso[2]=pfCandNeutralIso03Ele[theDenom1];
    
    chaPfIso[3]=pfCandChargedIso04Ele[theDenom1];
    phoPfIso[3]=pfCandPhotonIso04Ele[theDenom1];
    neuPfIso[3]=pfCandNeutralIso04Ele[theDenom1];
    
    chaPfIso[4]=pfCandChargedIso05Ele[theDenom1];
    phoPfIso[4]=pfCandPhotonIso05Ele[theDenom1];
    neuPfIso[4]=pfCandNeutralIso05Ele[theDenom1];
    
    chaPfIso[5]=pfCandChargedIso06Ele[theDenom1];
    phoPfIso[5]=pfCandPhotonIso06Ele[theDenom1];
    neuPfIso[5]=pfCandNeutralIso06Ele[theDenom1];
    
    chaPfIso[6]=pfCandChargedIso07Ele[theDenom1];
    phoPfIso[6]=pfCandPhotonIso07Ele[theDenom1];
    neuPfIso[6]=pfCandNeutralIso07Ele[theDenom1];
    
    chaPfIso[7]=pfCandChargedDirIso04Ele[theDenom1];
    phoPfIso[7]=pfCandPhotonDirIso04Ele[theDenom1];
    neuPfIso[7]=pfCandNeutralDirIso04Ele[theDenom1];

    float trkIso[2], ecalIso[2], hcalIso[2];
    trkIso[0]=dr03TkSumPtEle[theDenom1];
    trkIso[1]=dr04TkSumPtEle[theDenom1];
    ecalIso[0]=dr03EcalRecHitSumEtEle[theDenom1];
    ecalIso[1]=dr04EcalRecHitSumEtEle[theDenom1];
    hcalIso[0]=dr03HcalTowerSumEtFullConeEle[theDenom1];
    hcalIso[1]=dr04HcalTowerSumEtFullConeEle[theDenom1];
    
    bool isEleEB= anaUtils.fiducialFlagECAL(fiducialFlagsEle[theDenom1], isEB);
    bool isEleEE= anaUtils.fiducialFlagECAL(fiducialFlagsEle[theDenom1], isEE);

    // fill the reduced tree
    myOutIDTree->fillVariables(eleopout,eopout,eop,HoE,deta,dphi,s9s25,s1s9,see,spp,fbrem,
			       nbrems,misshits,dcot,dist,pt,eta,charge,phiwidth,etawidth,
			       oneoveremoneoverp,eledeta,d0,ip3d,ip3ds,kfhits,kflayers,kfchi2,e1x5e5x5,ecalseed,matchConv,
                               isEleEB,isEleEE);
    myOutIDTree->fillVariables2(detacalo, dphicalo, sep, dz, gsfchi2, emax, etop, ebottom, eleft, eright,
                                e2nd, e2x5right, e2x5left, e2x5top, e2x5bottom, 
                                e2x5max, e1x5, e2x2, e3x3, e5x5, r9, nclu,
                                phi, scenergy, scrawenergy, scesenergy,eseedopin);    
    myOutIDTree->fillCluterInfos(EleSCEt,EleSCEta,EleSCPhi,EtaSeed,PhiSeed,ESeed,IEtaSeed,IPhiSeed,EtaCrySeed,PhiCrySeed,IEtaCrySeed,IPhiCrySeed);
    myOutIDTree->fillIsolations(trkIso,ecalIso,hcalIso,
                                pfCombinedIsoEle[theDenom1],
                                chaPfIso, neuPfIso, phoPfIso);
    myOutIDTree->fillFakeRateDenomBits(ptLeadingJet,isDenomFake_HwwEgamma(theDenom1),isDenomFake_smurfs(theDenom1));
    myOutIDTree->fillMore(nPV,rhoFastjet,hwwbdts,newhwwbdts,hzzbdts,pfmva,lh);
    myOutIDTree->fillTrackMomenta(pcomb,pmodegsf,pmeangsf,pmeankf,pterrorgsf,pterrorkf,perrorele);
    myOutIDTree->fillCiCBasedIDBits(cic);
    myOutIDTree->fillRunInfos(runNumber, lumiBlock, eventNumber, nPU, -1);
    myOutIDTree->store();
    
  } // loop events

  // saving the counters
  sprintf(filename,"%sCounters.root",outname);
  myCounter.Save(filename,"recreate");
  
  // saving the output tree
  myOutKineTree -> save();
  myOutIDTree   -> save();

}

// denominator for fake rate: for HtoWW, egamma triggers
int FakeElectronSelector::isDenomFake_HwwEgamma(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);

  // match with the HLT firing candidates
  bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
  if (!HLTmatch) isGoodDenom = false;
  
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
int FakeElectronSelector::isDenomFake_smurfs(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // acceptance for the fake electron
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )        isGoodDenom = false;

  // match with the HLT firing candidates
  if(_isData) {
    bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
    if (!HLTmatch) isGoodDenom = false;
  }

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
