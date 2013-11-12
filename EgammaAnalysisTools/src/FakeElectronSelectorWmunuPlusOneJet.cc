
#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/FakeElectronSelectorWmunuPlusOneJet.hh"
#include "CommonTools/include/Counters.hh"
#include "EgammaAnalysisTools/include/SimpleCutsIDSelector.hh"
#include "EgammaAnalysisTools/include/eIDCiChzzSelector.hh"

using namespace bits;
using namespace std;

FakeElectronSelectorWmunuPlusOneJet::FakeElectronSelectorWmunuPlusOneJet(TTree *tree)
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
  myCounter.AddVar("wmuon");
  myCounter.AddVar("hltmatch");
  myCounter.AddVar("zveto");
  myCounter.AddVar("wmt");
  myCounter.AddVar("onejet");
  myCounter.AddVar("bveto");
  myCounter.AddVar("denom");
}

FakeElectronSelectorWmunuPlusOneJet::~FakeElectronSelectorWmunuPlusOneJet() { }


void FakeElectronSelectorWmunuPlusOneJet::Loop(const char *outname) {

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
  cout << "using muon triggers" << endl;
  requiredTriggers.push_back("HLT_Mu15_v2");
  requiredTriggers.push_back("HLT_Mu24_v");
  requiredTriggers.push_back("HLT_Mu30_v");
  requiredTriggers.push_back("HLT_Mu40_v");
  requiredTriggers.push_back("HLT_Mu40_eta2p1_v");
  requiredTriggers.push_back("HLT_IsoMu17_v");
  requiredTriggers.push_back("HLT_IsoMu20_v");
  requiredTriggers.push_back("HLT_IsoMu24_v");
  requiredTriggers.push_back("HLT_IsoMu24_eta2p1_v");

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
    // if ( _isData && !passedHLT ) continue;
    myCounter.IncrVar("trigger",1);

    // electrons passing the denominator selection (removed to get it unbiased. Set a bit at the end)
    SimpleCutsIDSelector aCutSel;
    Utils anaUtils;
    std::vector<int> tightMuons, looseMuons; 
    for(int imu=0; imu<nMuon; imu++) {

      float eta = etaMuon[imu];
      float pt = GetPt(pxMuon[imu],pyMuon[imu]);
      float tkiso = sumPt03Muon[imu]/pt;
      float ecalIsolAbs = 0.0;
      if ( fabs(eta)<1.479 ) ecalIsolAbs = max(0.0,emEt03Muon[imu]-1.0);
      else ecalIsolAbs = emEt03Muon[imu];
      float ecaliso = ecalIsolAbs/pt;
      float hcaliso = hadEt03Muon[imu]/pt;

      // use the same iolation as WP70 for electrons
      bool isGoodWMu = pt>25. && 
	aCutSel.outputiso(eta,tkiso,ecaliso,hcaliso,SimpleCutsIDSelector::kWP70);
      
      bool isLooseMu = aCutSel.outputiso(eta,tkiso,ecaliso,hcaliso,SimpleCutsIDSelector::kWP95);

      if (isGoodWMu) tightMuons.push_back(imu);
      if (isLooseMu) looseMuons.push_back(imu);
    }

    // event selection: at least one candidate for W->munu
    if (tightMuons.size()==0) continue;
    myCounter.IncrVar("wmuon",1);    

    // best muon as W candidate
    std::pair<int,int> bestpair = getBestGoodMuonPair(tightMuons);
    int theWMuon(bestpair.first);
    TLorentzVector wMuonP4;
    wMuonP4.SetXYZT(pxMuon[theWMuon],pyMuon[theWMuon],pzMuon[theWMuon],energyMuon[theWMuon]);

    // match with the HLT firing candidates
    if(_isData) {
      bool HLTmatch = triggerMatch(wMuonP4.Eta(),wMuonP4.Phi(),0.2);
      // if (!HLTmatch) continue;
    }
    myCounter.IncrVar("hltmatch",1);    

    // veto the Z on WP70 x WP95
    bool zmass = false;
    for(int iloose=0;iloose<(int)looseMuons.size() && iloose!=theWMuon && !zmass; ++iloose) {
      TLorentzVector looseP4;
      looseP4.SetXYZT(pxMuon[iloose],pyMuon[iloose],pzMuon[iloose],energyMuon[iloose]);
      float mass = (wMuonP4+looseP4).M();
      if (mass>60. && mass<120.) zmass = true;  
    }
    
    if(zmass) continue;
    myCounter.IncrVar("zveto",1);    

    if (theWMuon==-1) cout << "sanity check: impossibile!" << endl;

    // MT cut
    TVector3 p3Met(pxPFMet[0],pyPFMet[0],0.0);
    float WmT = sqrt(2*wMuonP4.Pt()*p3Met.Pt()*(1-cos(wMuonP4.Vect().Angle(p3Met))) );      
    if(WmT<50) continue;
    myCounter.IncrVar("wmt",1);  


    float maxTCHE = -1000.;
    float ptLeadingJet = -1;
    std::vector<int> jets30;
    for ( int jet=0; jet<nAK5PFPUcorrJet; jet++ ) {
      TVector3 p3Jet(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet]);
      if ( fabs(p3Jet.DeltaPhi(wMuonP4.Vect())) < 0.3 ) continue;
      if ( p3Jet.Pt()>30 && fabs(p3Jet.Eta())<3.0 ) jets30.push_back(jet); 
      if ( p3Jet.Pt()>15 && fabs(p3Jet.Eta())<3.0 && trackCountingHighEffBJetTagsAK5Jet[jet]>maxTCHE ) maxTCHE= trackCountingHighEffBJetTagsAK5Jet[jet];
      if ( p3Jet.Pt()>ptLeadingJet ) ptLeadingJet = p3Jet.Pt();
    }

    // require <=1 jet with pt > 30 GeV
    if(jets30.size()>1) continue;
    myCounter.IncrVar("onejet",1);

    // veto on max btag in jets with pt > 15 GeV
    if(maxTCHE>2.1) continue;
    myCounter.IncrVar("bveto",1);  

    // look for probes among the reco electrons
    std::vector<int> probes;
    for(int iele=0; iele<nEle; iele++) {
      TVector3 probeP3;
      probeP3.SetXYZ(pxEle[iele],pyEle[iele],pzEle[iele]);
      // exclude possible brem candidates in a cone 0.3 from W muon
      if(fabs(probeP3.Eta())<2.5 && probeP3.Pt()>5. && fabs(probeP3.DeltaPhi(wMuonP4.Vect()))>0.3) probes.push_back(iele);
    }

    // take the highest pt probe
    std::pair<int,int> bestprobes = getBestGoodElePair(probes);     
    int probe(bestprobes.first);
    TLorentzVector probeP4;
    probeP4.SetXYZT(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);

    // at least one probe candidate
    if(probe<0) continue;
    myCounter.IncrVar("denom",1);
    
    float mass = -1.;

    // fill the kine tree - after HLT and denominator
    myOutKineTree -> fill( p3Met.Pt(), WmT, mass, wMuonP4.Pt(), probeP4.Pt(), fabs(wMuonP4.Vect().DeltaPhi(probeP4.Vect())) );
    myOutKineTree -> store();
    
    // fill the denominator: take only the highest pT denominator candidate
    float etaFake = fabs(probeP4.Eta()); 
    float etFake  = probeP4.Pt();

    // some eleID variables
    float HoE, s1s9, s9s25, phiwidth, etawidth, deta, dphi, fbrem, see, spp, eleopout, eopout, eop, eseedopin, nbrems, recoFlag, EleSCEta, EleSCPhi;
    float EleSCEt,EtaSeed,PhiSeed,ESeed,IEtaSeed,IPhiSeed,EtaCrySeed,PhiCrySeed,IEtaCrySeed,IPhiCrySeed;
    float oneoveremoneoverp, eledeta, d0, ip3d, ip3ds, kfhits, kflayers, kfchi2, e1x5e5x5, dcot, dist;
    float detacalo, dphicalo, sep, dz, gsfchi2, emax, etop, ebottom, eleft, eright,
      e2nd, e2x5right, e2x5left, e2x5top, e2x5bottom, 
      e2x5max, e1x5, e2x2, e3x3, e5x5, r9, nclu,
      scenergy, scrawenergy, scesenergy;

    int gsfTrack = gsfTrackIndexEle[probe];
    int kfTrack = trackIndexEle[probe];

    // different p estimations
    float pcomb=probeP4.Vect().Mag();
    TVector3 p3ModeGsf(pxModeGsfTrack[gsfTrack],pyModeGsfTrack[gsfTrack],pzModeGsfTrack[gsfTrack]);
    float pmodegsf=p3ModeGsf.Mag();
    TVector3 p3MeanGsf(pxGsfTrack[gsfTrack],pyGsfTrack[gsfTrack],pzGsfTrack[gsfTrack]);
    float pmeangsf=p3MeanGsf.Mag();
    TVector3 p3MeanKf(pxTrack[kfTrack],pyTrack[kfTrack],pzTrack[kfTrack]);
    float pmeankf=p3MeanKf.Mag();
    float pterrorgsf = ptErrorGsfTrack[gsfTrack];
    float pterrorkf = (kfTrack>-1) ? ptErrorTrack[kfTrack] : -1.0;
    float perrorele = trackMomentumErrorEle[probe];

    double gsfsign   = (-eleDxyPV(probe,0) >=0 ) ? 1. : -1.;
    bool matchConv = hasMatchedConversionEle[probe];

    d0 = gsfsign * transvImpactParGsfTrack[gsfTrack];
    dz = eleDzPV(probe,0);
    ip3d = gsfsign * impactPar3DGsfTrack[gsfTrack];
    ip3ds = ip3d/impactPar3DErrorGsfTrack[gsfTrack];
    kfchi2 = (kfTrack>-1) ? trackNormalizedChi2Track[kfTrack] : 0.0;
    kflayers = (kfTrack>-1) ? trackerLayersWithMeasurementTrack[kfTrack] : -1.0;
    kfhits = (kfTrack>-1) ? trackValidHitsTrack[kfTrack] : -1.0;
    gsfchi2 = trackNormalizedChi2GsfTrack[gsfTrack];
    int misshits = expInnerLayersGsfTrack[gsfTrack];
    dcot = convDistEle[probe];
    dist = convDcotEle[probe];
    bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[probe], isEcalDriven);
    int ecalseed;
    HoE = hOverEEle[probe];
    eledeta = deltaEtaEleClusterTrackAtCaloEle[probe];
    deta = deltaEtaAtVtxEle[probe];
    dphi = deltaPhiAtVtxEle[probe];
    detacalo = deltaEtaAtCaloEle[probe];
    dphicalo = deltaPhiAtCaloEle[probe];
    fbrem = fbremEle[probe];
    nbrems = nbremsEle[probe];
    eleopout = eEleClusterOverPoutEle[probe];
    eopout = eSeedOverPoutEle[probe];
    eop = eSuperClusterOverPEle[probe];
    if(ecaldriven) {
      ecalseed = 1;
      int sc = superClusterIndexEle[probe];
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
      oneoveremoneoverp = 1./energySC[sc]  - 1./probeP4.Vect().Mag();
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
      int sc = PFsuperClusterIndexEle[probe];
      if(sc>-1) {
	float seedEnergy = seedClusterEnergyPFSC[sc];
        nclu = float(nBCPFSC[sc]);
        s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
        s1s9 = eMaxPFSC[sc]/eMaxPFSC[sc];
        e1x5e5x5 = (e5x5PFSC[sc] - e1x5PFSC[sc])/e5x5PFSC[sc];
        phiwidth = phiWidthPFSC[sc];
        etawidth = etaWidthPFSC[sc];
        see = sqrt(covIEtaIEtaPFSC[sc]);
        sep = covIEtaIPhiPFSC[sc]/(sqrt(covIEtaIEtaPFSC[sc])*sqrt(covIPhiIPhiPFSC[sc]));
        spp = sqrt(covIPhiIPhiPFSC[sc]);
        oneoveremoneoverp = 1./energyPFSC[sc]  - 1./probeP4.Vect().Mag();
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

    float pt = probeP4.Pt();
    float eta = probeP4.Eta();
    float phi = probeP4.Phi();
    int charge = chargeEle[probe];

    // CiC...
    eIDCiChzzSelector cicsel;
    int cic[5];
    int iecal = (fabs(eta)<1.479) ? 0 : 1;
    for(int i=0; i<5; i++) cic[i] = cicsel.ElectronId_V03(pt,EleSCEta,see,eop,eseedopin,fbrem,
							  dr03TkSumPtEle[probe],
							  dr03EcalRecHitSumEtEle[probe] - rhoFastjet * Aeff_ecal_dr03[iecal],
							  dr03HcalTowerSumEtEle[probe] - rhoFastjet * Aeff_hcal_dr03[iecal],
							  d0,misshits,deta,dphi,HoE,dcot,
							  dist,!ecalseed,i,false);
    
    // some MVAs...
    float pfmva = pflowMVAEle[probe];
    float lh=eleIdLikelihoodEle[probe];
    float hwwbdts[2];
    hwwbdts[0] = eleBDT(fMVAHWW,probe);
    hwwbdts[1] = eleBDT(fMVAHWWWithIso,probe);
    float hzzbdts[4];
    hzzbdts[0] = eleBDT(fMVAHZZDanV0,probe);
    hzzbdts[1] = eleBDT(fMVAHZZSiV0,probe);
    hzzbdts[2] = eleBDT(fMVAHZZSiV1,probe);
    // hzzbdts[3] = eleBDT(fMVAHZZSiDanV2,probe);
    hzzbdts[3] = mvaidnontrigEle[probe];
    float newhwwbdts[4];
    newhwwbdts[0] = eleBDT(fMVAHWWDanV0,probe);
    newhwwbdts[1] = eleBDT(fMVAHWWSiV0,probe);
    newhwwbdts[2] = eleBDT(fMVAHWWSiV1,probe);
    // newhwwbdts[3] = eleBDT(fMVAHWWSiDanV2,probe);
    newhwwbdts[3] = mvaidtrigEle[probe];

    // isolations
    float chaPfIso[8], phoPfIso[8], neuPfIso[8];
    chaPfIso[0]=pfCandChargedIso01Ele[probe];
    phoPfIso[0]=pfCandPhotonIso01Ele[probe];
    neuPfIso[0]=pfCandNeutralIso01Ele[probe];
    
    chaPfIso[1]=pfCandChargedIso02Ele[probe];
    phoPfIso[1]=pfCandPhotonIso02Ele[probe];
    neuPfIso[1]=pfCandNeutralIso02Ele[probe];
    
    chaPfIso[2]=pfCandChargedIso03Ele[probe];
    phoPfIso[2]=pfCandPhotonIso03Ele[probe];
    neuPfIso[2]=pfCandNeutralIso03Ele[probe];
    
    chaPfIso[3]=pfCandChargedIso04Ele[probe];
    phoPfIso[3]=pfCandPhotonIso04Ele[probe];
    neuPfIso[3]=pfCandNeutralIso04Ele[probe];
    
    chaPfIso[4]=pfCandChargedIso05Ele[probe];
    phoPfIso[4]=pfCandPhotonIso05Ele[probe];
    neuPfIso[4]=pfCandNeutralIso05Ele[probe];
    
    chaPfIso[5]=pfCandChargedIso06Ele[probe];
    phoPfIso[5]=pfCandPhotonIso06Ele[probe];
    neuPfIso[5]=pfCandNeutralIso06Ele[probe];
    
    chaPfIso[6]=pfCandChargedIso07Ele[probe];
    phoPfIso[6]=pfCandPhotonIso07Ele[probe];
    neuPfIso[6]=pfCandNeutralIso07Ele[probe];
    
    chaPfIso[7]=pfCandChargedDirIso04Ele[probe];
    phoPfIso[7]=pfCandPhotonDirIso04Ele[probe];
    neuPfIso[7]=pfCandNeutralDirIso04Ele[probe];

    float trkIso[2], ecalIso[2], hcalIso[2];
    trkIso[0]=dr03TkSumPtEle[probe];
    trkIso[1]=dr04TkSumPtEle[probe];
    ecalIso[0]=dr03EcalRecHitSumEtEle[probe];
    ecalIso[1]=dr04EcalRecHitSumEtEle[probe];
    hcalIso[0]=dr03HcalTowerSumEtFullConeEle[probe];
    hcalIso[1]=dr04HcalTowerSumEtFullConeEle[probe];
    
    bool isEleEB= anaUtils.fiducialFlagECAL(fiducialFlagsEle[probe], isEB);
    bool isEleEE= anaUtils.fiducialFlagECAL(fiducialFlagsEle[probe], isEE);

    // fill the reduced tree
    myOutIDTree->fillVariables(eleopout,eopout,eop,HoE,deta,dphi,s9s25,s1s9,see,spp,fbrem,
			       nbrems,misshits,dcot,dist,pt,eta,charge,phiwidth,etawidth,
			       oneoveremoneoverp,eledeta,d0,ip3d,ip3ds,kfhits,kflayers,kfchi2,e1x5e5x5,ecalseed,matchConv,
                               isEleEB,isEleEE);
    myOutIDTree->fillVariables2(detacalo, dphicalo, sep, dz, gsfchi2, emax, etop, ebottom, eleft, eright,
                                e2nd, e2x5right, e2x5left, e2x5top, e2x5bottom, 
                                e2x5max, e1x5, e2x2, e3x3, e5x5, r9, nclu,
                                phi, scenergy, scrawenergy, scesenergy, eseedopin);    
    myOutIDTree->fillCluterInfos(EleSCEt,EleSCEta,EleSCPhi,EtaSeed,PhiSeed,ESeed,IEtaSeed,IPhiSeed,EtaCrySeed,PhiCrySeed,IEtaCrySeed,IPhiCrySeed);
    myOutIDTree->fillIsolations(trkIso,ecalIso,hcalIso,
                                pfCombinedIsoEle[probe],
                                chaPfIso, neuPfIso, phoPfIso);
    myOutIDTree->fillFakeRateDenomBits(ptLeadingJet,isDenomFake_HwwEgamma(probe),isDenomFake_smurfs(probe));
    myOutIDTree->fillMore(nPV,rhoFastjet,hwwbdts,newhwwbdts,hzzbdts,pfmva,lh);
    myOutIDTree->fillTrackMomenta(pcomb,pmodegsf,pmeangsf,pmeankf,pterrorgsf,pterrorkf,perrorele);
    myOutIDTree->fillCiCBasedIDBits(cic);
    myOutIDTree->fillRunInfos(runNumber, lumiBlock, eventNumber, nPU, -1);
    myOutIDTree->store();
    
  } // loop events

  // saving the counters
  sprintf(filename,"%sCounters.root",outname);
  myCounter.Save(filename,"recreate");
  myCounter.Draw();
  
  // saving the output tree
  myOutKineTree -> save();
  myOutIDTree   -> save();

}

// two highest pT muons
std::pair<int,int> FakeElectronSelectorWmunuPlusOneJet::getBestGoodMuonPair(std::vector<int> goodMuons) {
  
  int theMuon1=-1;
  int theMuon2=-1;
  float maxPt1=-1000.;
  float maxPt2=-1001.;
  for(int iMuon=0;iMuon<goodMuons.size();iMuon++) {
    int muonIndex = goodMuons[iMuon];
    TVector3 pMuon(pxMuon[muonIndex],pyMuon[muonIndex],pzMuon[muonIndex]);
    float thisPt=pMuon.Pt();
    if (thisPt>maxPt1 && thisPt>maxPt2){ maxPt2 = maxPt1; maxPt1 = thisPt; theMuon2 = theMuon1; theMuon1 = muonIndex; }
    if (thisPt<maxPt1 && thisPt>maxPt2){ maxPt2 = thisPt; theMuon2 = muonIndex; }
  }
  return std::make_pair(theMuon1,theMuon2);
}
