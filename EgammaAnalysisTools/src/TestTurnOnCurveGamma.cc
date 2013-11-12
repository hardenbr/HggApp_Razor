#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/TestTurnOnCurveGamma.hh"
#include "CommonTools/include/Counters.hh"
#include <CommonTools/include/TriggerMask.hh>
#include "EgammaAnalysisTools/include/TurnOnTreeGamma.hh"

using namespace bits;
using namespace std;

TestTurnOnCurveGamma::TestTurnOnCurveGamma(TTree *tree)
  : Egamma(tree) {
  
  _isData = true;    
  
  // to read good run list
  if (_isData) {
    std::string goodRunGiasoneFile = "config/json/ABCD_11DecExcluded.json"; 
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }

  // single electron efficiency - simple cuts 2011 (ours and smurfs) - chiara, vedi cosa usare
  EgammaCutBasedIDLowPtWPs.push_back("WP95");       // 0
  EgammaCutBasedIDLowPtWPs.push_back("WP90");       // 1
  EgammaCutBasedIDLowPtWPs.push_back("WP85");       // 2
  EgammaCutBasedIDLowPtWPs.push_back("WP80");       // 3
  EgammaCutBasedIDLowPtWPs.push_back("WP70");       // 4
  EgammaCutBasedIDLowPtWPs.push_back("WP70Smurf");  // 5

  // single electron efficiency
  for (int i=0;i<EgammaCutBasedIDLowPtWPs.size();++i) {
    CutBasedEleIDSelector aSelector;
    char configDir[50];
    sprintf(configDir,"config/%s",EgammaCutBasedIDLowPtWPs[i].c_str());
    std::cout << "===== Configuring " <<  EgammaCutBasedIDLowPtWPs[i] << " ElectronID ==========" << std::endl;
    aSelector.ConfigureNoClass(configDir);
    EgammaCutBasedIDLowPt.push_back(aSelector);
  }
}

TestTurnOnCurveGamma::~TestTurnOnCurveGamma() { }

void TestTurnOnCurveGamma::measureTurnOnGamma(const char *outname) {

  // output tree
  char treename[200];
  sprintf(treename,"%s-tree.root",outname);
  TurnOnTreeGamma reducedTree(treename);
  reducedTree.addTandPinfo();
  reducedTree.addHLTmatchInfo();
  reducedTree.addMore();
  reducedTree.addRunInfos();

  // trigger: here is the lowest pT HLT path 
  std::vector<std::string> mask;
  // mask.push_back("HLT_Photon20_CaloIdVL_IsoL_v");      
  // mask.push_back("HLT_Photon30_CaloIdVL_IsoL_v");      
  // mask.push_back("HLT_Photon50_CaloIdVL_IsoL_v");      
  mask.push_back("HLT_Photon75_CaloIdVL_IsoL_v");      

  // json
  unsigned int lastLumi = 0;
  unsigned int lastRun  = 0;

  // loop on events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // reload trigger mask (the bit)
    setRequiredTriggers(mask);
    reloadTriggerMask(true);

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
    
    // only events passing the lowest pT threshold single photon trigger path
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if (!passedHLT ) continue;   

    // check if the event is a Z one - to have a pure enough sample (can be used offline)
    bool isAZ       = true;
    bool isAZ_tight = true;
    
    // a) at least 2 reco electrons
    if (nEle<2) { 
      isAZ = false; 
      isAZ_tight = false; 
    }

    // b) best tag-probe pair = mee closest to Z mass
    float minpull = 100000.;
    float mass    = 1000.;
    float okmass  = 1000.;
    int theEle1   = -1;
    int theEle2   = -1;
    float ptEle1  = 999.;
    float ptEle2  = 999.;
    float etaEle1 = 999.;
    float etaEle2 = 999.;

    for(int iele1=0; iele1<nEle; iele1++) {
      TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);
      
      for(int iele2=iele1+1; iele2<nEle; iele2++) {
        TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);

	// ECAL driven electrons only
	int sc1 = superClusterIndexEle[iele1];
	int sc2 = superClusterIndexEle[iele2];
	if ( sc1 < 0 || sc2 < 0 ) continue;

	// eta cuts
	if ( fabs(etaSC[sc1])>1.4442 && fabs(etaSC[sc1])<1.566 ) continue;  
	if ( fabs(etaSC[sc2])>1.4442 && fabs(etaSC[sc2])<1.566 ) continue;  
	if ( fabs(etaSC[sc1])>2.5 ) continue;  
	if ( fabs(etaSC[sc2])>2.5 ) continue;  

	// pT cuts
        if ( electron1.Pt()<20. ) continue;
        if ( electron2.Pt()<20. ) continue;

	// best invariant mass
        mass = (electron1+electron2).M();
        float pull=fabs(mass-91.1876);
	if(pull < minpull) {
	  okmass  = mass;
          minpull = pull;
          theEle1 = iele1;     
          theEle2 = iele2;     
	  ptEle1  = electron1.Pt();
	  ptEle2  = electron2.Pt();
	  etaEle1 = etaSC[sc1];
	  etaEle2 = etaSC[sc2];
        }
      }
    }
    if ( okmass>120 || okmass<60 ) { 
      isAZ = false; 
      isAZ_tight = false; 
    }

    // c) one of the two electrons should be tight enough (a sort of tag)
    bool tagIdentified1, tagIsolated1, tagConvRej1;
    bool tagIdentified2, tagIsolated2, tagConvRej2;
    tagIdentified1 = tagIsolated1 = tagConvRej1 = false;
    tagIdentified2 = tagIsolated2 = tagConvRej2 = false;
    if ( isAZ ) {
      isEleID(&EgammaCutBasedIDLowPt[5],theEle1,&tagIdentified1,&tagIsolated1,&tagConvRej1);
      isEleID(&EgammaCutBasedIDLowPt[5],theEle2,&tagIdentified2,&tagIsolated2,&tagConvRej1);
      if ( (!tagIdentified1 || !tagIsolated1) && (!tagIdentified2 || !tagIsolated2) ) isAZ_tight = false;
    }
    bool tag1 = tagIdentified1 && tagIsolated1;
    bool tag2 = tagIdentified2 && tagIsolated2;
    
    // if wanted for checks - chiara: per prove lo metto nel tree
    // if ( !isAZ )       continue;
    // if ( !isAZ_tight ) continue;


    // at least 1 reconstructed photon candidate passing the full selection and matching the HLT object
    int maxPtGoodGamma = -1;
    float ptGoodGamma = -999.;
    for(int ipho=0; ipho<nPho; ipho++) {

      // eta cuts
      bool isEB = false;
      if ( fabs(etaPho[ipho])>1.4442 && fabs(etaPho[ipho])<1.566 ) continue;  
      if ( fabs(etaPho[ipho])>2.5 ) continue;  
      if ( fabs(etaPho[ipho])<1.479) isEB = true;
      
      // pT cuts
      TLorentzVector photonP4(pxPho[ipho],pyPho[ipho],pzPho[ipho],energyPho[ipho]);
      if ( photonP4.Pt()<20. ) continue;
      
      // match with HLT candidate
      // bool matchGamma = triggerMatch(etaPho[ipho],phiPho[ipho],0.2);
      // if (!matchGamma) continue;                   // chiara: per prove lo metto nel tree
      
      // pass the full offline selection
      // bool isFullSelGamma = isPhotonID(ipho);      // chiara: per prove lo metto nel tree
      // if (!isFullSelGamma) continue;
      
      if ( photonP4.Pt() > ptGoodGamma ){
	ptGoodGamma = photonP4.Pt();
	maxPtGoodGamma = ipho;
      }
    }
	
    // one gamma is enough
    if (maxPtGoodGamma<0) continue;
    bool isGoodEB = false;
    if (fabs(etaPho[maxPtGoodGamma])<1.479) isGoodEB = true; 
    // for the tree
    int matchGoodGamma = triggerMatch(etaPho[maxPtGoodGamma],phiPho[maxPtGoodGamma],0.2);
    int goodGoodGamma  = isPhotonID(maxPtGoodGamma);

    // consider the HLT object firing the N-1 HLT path
    bool hltMatchAboveThreshold_20 = triggerMatchThreshold(20);
    bool hltMatchAboveThreshold_30 = triggerMatchThreshold(30);
    bool hltMatchAboveThreshold_50 = triggerMatchThreshold(50);
    bool hltMatchAboveThreshold_75 = triggerMatchThreshold(75);
    bool hltMatchAboveThreshold_90 = triggerMatchThreshold(90);
    
    // fill the reduced tree 
    reducedTree.fillVariables(ptGoodGamma, etaPho[maxPtGoodGamma], phiPho[maxPtGoodGamma], goodGoodGamma, matchGoodGamma);
    reducedTree.fillTandPinfo(okmass, ptEle1, ptEle2, etaEle1, etaEle2, tag1, tag2, isAZ, isAZ_tight);
    reducedTree.fillHLTmatchInfo(hltMatchAboveThreshold_20, hltMatchAboveThreshold_30, hltMatchAboveThreshold_50, hltMatchAboveThreshold_75, hltMatchAboveThreshold_90);
    reducedTree.fillMore(nPV,rhoFastjet);
    reducedTree.fillRunInfos(runNumber, lumiBlock, eventNumber);

    // store
    reducedTree.store();

  } // loop on events

  reducedTree.save();
}

// apply the full cut based photonID selection (medium WP)                                                                                         
// https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012                                                                                       
bool TestTurnOnCurveGamma::isPhotonID(int theGamma) {

  bool isFullySel = true;

  // supercluster associated to the photon
  int theSc = superClusterIndexPho[theGamma];
  if ( theSc<0 ) return false;
  float etasc = etaSC[theSc];
  float see   = sqrt(covIEtaIEtaSC[theSc]);

  // photon 4vector
  TLorentzVector photonP4(pxPho[theGamma],pyPho[theGamma],pzPho[theGamma],energyPho[theGamma]);
  float pt = photonP4.Pt();

  // EA corrections
  float EA_chargedH[7] = { 0.012, 0.010, 0.014, 0.012, 0.016, 0.020, 0.012};
  float EA_neutralH[7] = { 0.030, 0.057, 0.039, 0.015, 0.024, 0.039, 0.072};
  float EA_photons[7]  = { 0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};
  
  // for effective area calculation       
  int theEAregion = effectiveAreaRegion(etasc);
  if (theEAregion>6) return false;
  
  // corrected isolations - chiara: va fatto cosi', su ciascuno? rhoJetsFastJet e' quello ok? cono 03?  
  // chiara: problema, credo che il cono sia 0.4
  float rhoCorrCharged = chargedHadronIsoPho[theGamma] - rhoJetsFastJet*EA_chargedH[theEAregion];                                               
  float rhoCorrNeutral = neutralHadronIsoPho[theGamma] - rhoJetsFastJet*EA_neutralH[theEAregion];                                              
  float rhoCorrPhoton  = photonIsoPho[theGamma] - rhoJetsFastJet*EA_photons[theEAregion];                                                
  if (rhoCorrCharged<0) rhoCorrCharged = 0.;                                                                                                       
  if (rhoCorrNeutral<0) rhoCorrNeutral = 0.;                                                                                                       
  if (rhoCorrPhoton<0)  rhoCorrPhoton  = 0.;  
  
  if(hOverEPho[theGamma]>0.05)  isFullySel = false;                                                                                    
  if (theEAregion<2) {  // EB                                                                                                                      
    if (see>0.011) isFullySel = false;                                                                                    
    if (rhoCorrCharged > 1.5) isFullySel = false;                                                                                    
    if (rhoCorrNeutral > 1.0 + 0.04*pt)  isFullySel = false;                                                                  
    if (rhoCorrPhoton  > 0.7 + 0.005*pt) isFullySel = false;                                                                  
  } else {     // EE                                                                                                                               
    if (see>0.033) isFullySel = false;                                                                                    
    if (rhoCorrCharged > 1.2) isFullySel = false;                                                                                    
    if (rhoCorrNeutral > 1.5 + 0.04*pt)  isFullySel = false;                                                                  
    if (rhoCorrPhoton  > 1.0 + 0.005*pt) isFullySel = false;                                                                  
  }
  
  return isFullySel;
}

// for effective area calculation                                                                                                                      
int TestTurnOnCurveGamma::effectiveAreaRegion(float theEta) {
  
  int theEAregion = 999;
  if (fabs(theEta)<1.) theEAregion = 0;
  if (fabs(theEta)<1.479 && fabs(theEta)>1.)    theEAregion = 1;
  if (fabs(theEta)<2.    && fabs(theEta)>1.479) theEAregion = 2;
  if (fabs(theEta)<2.2   && fabs(theEta)>2.0)   theEAregion = 3;
  if (fabs(theEta)<2.3   && fabs(theEta)>2.2)   theEAregion = 4;
  if (fabs(theEta)<2.4   && fabs(theEta)>2.3)   theEAregion = 5;
  if (fabs(theEta)>2.4) theEAregion = 6;
  return theEAregion;
}

void TestTurnOnCurveGamma::isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;
  
  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  
  TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  float e1, e4SwissCross, fidFlagSC, seedRecHitFlag, seedTime, seedChi2;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];
  deta = deltaEtaAtVtxEle[eleIndex];
  dphiin = deltaPhiAtVtxEle[eleIndex];

  dphiout = deltaPhiAtCaloEle[eleIndex];
  fbrem = fbremEle[eleIndex];
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
    e1 = eMaxSC[sc];
    e4SwissCross = e4SwissCrossSC[sc];
    fidFlagSC = fiducialFlagsEle[eleIndex];
    seedRecHitFlag = recoFlagSC[sc];
    seedTime = timeSC[sc];
    seedChi2 = chi2SC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxPFSC[sc];
      e4SwissCross = e4SwissCrossPFSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagPFSC[sc];
      seedTime = timePFSC[sc];
      seedChi2 = chi2PFSC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  float lh = eleIdLikelihoodEle[eleIndex];
  bool isEleEB= anaUtils.fiducialFlagECAL(fiducialFlagsEle[eleIndex], isEB);
  selector->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  selector->SetRecoFlag(recoFlagsEle[eleIndex]);
  selector->applyElectronIDOnPFlowElectrons(true);
  selector->SetHOverE( HoE );
  selector->SetS9S25( s9s25 );
  selector->SetNBrem( nbremsEle[eleIndex] );
  selector->SetDEta( deta );
  selector->SetDPhiIn( dphiin );
  selector->SetDPhiOut( dphiout );
  selector->SetBremFraction( fbrem );
  selector->SetSigmaEtaEta( see );
  selector->SetSigmaPhiPhi( spp );
  selector->SetEOverPout( eopout );
  selector->SetEOverPin( eop );
  selector->SetElectronClass ( classificationEle[eleIndex] );
  selector->SetLikelihood( lh );
  selector->SetEcalIsolation( (dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  selector->SetTrkIsolation( (dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  selector->SetHcalIsolation( (dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  float combinedIso = 0.0;
  if (isEleEB) combinedIso = dr03TkSumPtEle[eleIndex] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else combinedIso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  selector->SetCombinedIsolation( (combinedIso - rhoFastjet*TMath::Pi()*0.3*0.3) / pEle.Pt() );

  selector->SetCombinedPFIsolation( (pfCombinedIsoEle[eleIndex]) / pEle.Pt() );

  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );
  selector->SetHasMatchedConversion ( hasMatchedConversionEle[eleIndex] );

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputNoClassEleId();
  *isolOutput = selector->outputNoClassIso();
  *convRejOutput = selector->outputNoClassConv();
}
