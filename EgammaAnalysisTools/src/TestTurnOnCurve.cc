#include <iostream>
#include <math.h>

#include "TVector3.h"
#include "TH1F.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/TestTurnOnCurve.hh"
#include "CommonTools/include/Counters.hh"

using namespace bits;
using namespace std;

TestTurnOnCurve::TestTurnOnCurve(TTree *tree)
  : Egamma(tree) {
  
  _isData = true;    
  
  // to read good run list
  if (_isData) {
    std::string goodRunGiasoneFile = "/afs/cern.ch/user/c/crovelli/scratch0/Vecbos2010/HiggsAnalysisTools/config/json/goodCollisions2011.json"; 
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }
}

TestTurnOnCurve::~TestTurnOnCurve() { }

void TestTurnOnCurve::measureTurnOn(const char *outname) {
  
  TH1F *RelIsoDenomEB    = new TH1F( "RelIsoDenomEB",    "denominator relative iso, EB",          75, -0.5, 1. );
  TH1F *RelIsoDenomEE    = new TH1F( "RelIsoDenomEE",    "denominator relative iso, EE",          75, -0.5, 1. );
  TH1F *RelIsoNumEB      = new TH1F( "RelIsoNumEB",      "numerator relative iso - no match, EB", 75, -0.5, 1. );
  TH1F *RelIsoNumEE      = new TH1F( "RelIsoNumEE",      "numerator relative iso - no match, EE", 75, -0.5, 1. );
  TH1F *RelIsoNumMatchEB = new TH1F( "RelIsoNumMatchEB", "denominator relative iso - match, EB",  75, -0.5, 1. );
  TH1F *RelIsoNumMatchEE = new TH1F( "RelIsoNumMatchEE", "denominator relative iso - match, EE",  75, -0.5, 1. );
 

  // trigger: here are those path I want to test
  // requiredTriggers.push_back("HLT_Ele17_CaloIdL_CaloIsoVL_v2");
  requiredTriggers.push_back("HLT_Ele8_CaloIdL_CaloIsoVL_v2");


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

    Utils anaUtils;
    
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
    

    // event selection: trigger => here are my "denominator" HLT paths ( those I use to select events )
    std::vector<int> thisTriggerMask;
    std::vector<std::string> theRequiredTriggers;
    std::vector<int> m_theRequiredTriggers;
    theRequiredTriggers.push_back("HLT_Ele8_v2");  
    for (std::vector< std::string >::const_iterator fIter=theRequiredTriggers.begin();fIter!=theRequiredTriggers.end();++fIter){
      for(unsigned int i=0; i<nameHLT->size(); i++){
	if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) ){
	  thisTriggerMask.push_back( indexHLT[i] ) ;
	  break;
	}
      }
    }
    m_theRequiredTriggers = thisTriggerMask;
    bool hasPassedDenomHLT = anaUtils.getTriggersOR(m_theRequiredTriggers, firedTrg);
    if ( !hasPassedDenomHLT ) continue;   
    
    // to clean up
    TVector3 p3Met(pxPFMet[0],pyPFMet[0],0.0);
    if( p3Met.Pt() > 20 ) continue;        
    
    // electrons passing the emulation of HLT cuts + our isolation requirement
    std::vector<int> denomElectrons;
    for(int iele=0; iele<nEle; iele++) {
      bool isGoodDenom = isDenomTurnOn(iele);
      if (!isGoodDenom) continue;
      denomElectrons.push_back(iele);
    }

    // considering electrons passing the baseline: check our isolation definition for them
    for (int iele=0; iele<denomElectrons.size(); iele++){
      
      int theEle = denomElectrons[iele];
      TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);

      float combinedIso;
      bool isEleEB = false;
      if ( fabs(etaEle[theEle]) < 1.479 ) isEleEB = true;
      if ( isEleEB) combinedIso = dr03TkSumPtEle[theEle] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[theEle]-1.0) + dr03HcalTowerSumEtFullConeEle[theEle];
      if (!isEleEB) combinedIso = dr03TkSumPtEle[theEle] + dr03EcalRecHitSumEtEle[theEle] + dr03HcalTowerSumEtFullConeEle[theEle];
      float corrCombinedIso = (combinedIso - rhoFastjet*TMath::Pi()*0.3*0.3) / p3Ele.Pt();       
      // float corrCombinedIso = (dr03EcalRecHitSumEtEle[theEle])/p3Ele.Pt();
      
      if ( isEleEB ) RelIsoDenomEB -> Fill( corrCombinedIso );
      if (!isEleEB ) RelIsoDenomEE -> Fill( corrCombinedIso );

      // now we require the event to pass the trigger we want to study
      reloadTriggerMask(true);
      bool passedHLT = hasPassedHLT();
      if ( passedHLT ) {
	
	if ( isEleEB ) RelIsoNumEB -> Fill( corrCombinedIso );
	if (!isEleEB ) RelIsoNumEE -> Fill( corrCombinedIso );
	
	// match with the HLT firing candidates
	bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
	if (HLTmatch) {
	  if ( isEleEB ) RelIsoNumMatchEB -> Fill( corrCombinedIso );
	  if (!isEleEB ) RelIsoNumMatchEE -> Fill( corrCombinedIso );
	}	
      }
    } // denominators    

  } // loop on events


  char filename[200];
  sprintf(filename,"%s-turnon.root",outname);
  TFile myFile(filename, "RECREATE");
  RelIsoDenomEB    -> Write();
  RelIsoDenomEE    -> Write();
  RelIsoNumEB      -> Write();
  RelIsoNumEE      -> Write();
  RelIsoNumMatchEB -> Write();
  RelIsoNumMatchEE -> Write();
}

// offline requirements for the turn on curve
bool TestTurnOnCurve::isDenomTurnOn(int theEle) {
  
  Utils anaUtils;
  
  bool isGoodDenom = true;

  // acceptance  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  // if( p3Ele.Pt() < 15. )        isGoodDenom = false;
  if( p3Ele.Pt() < 20. )        isGoodDenom = false;
  
  // H/E - to emulate the trigger
  bool isEleEB = false;
  if ( fabs(etaEle[theEle]) <  1.479 ) isEleEB = true;
  if ( isEleEB && hOverEEle[theEle]>0.15) isGoodDenom = false;
  if (!isEleEB && hOverEEle[theEle]>0.10) isGoodDenom = false;
  
  // sigmaIetaIeta - to emulate the trigger
  bool isBarrelSc;
  int sc = superClusterIndexEle[theEle];
  if ( sc < 0 ) isGoodDenom = false;
  if ( fabs(etaSC[sc]) <  1.479 ) isBarrelSc = true;
  if ( fabs(etaSC[sc]) >= 1.479 ) isBarrelSc = false;
  if ( isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.014 ) isGoodDenom = false;
  if (!isBarrelSc && sqrt(covIEtaIEtaSC[sc])>0.035 ) isGoodDenom = false;
      
  return isGoodDenom;
}
    
