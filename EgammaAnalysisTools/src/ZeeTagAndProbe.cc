#include <iostream>
#include <string> 
#include <math.h>

#include "TLorentzVector.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "EgammaAnalysisTools/include/ZeeTagAndProbe.hh"
#include "EgammaAnalysisTools/include/eIDCiChzzSelector.hh"

#include "cajun/json/reader.h"
#include "cajun/json/elements.h"
#include <CommonTools/include/TriggerMask.hh>

using namespace bits;
using namespace std;

ZeeTagAndProbe::ZeeTagAndProbe(TTree *tree)
  : Egamma(tree) {

  jsonFile = "";
  lastFile = "";
  
  std::string fileCuts("config/ZeeTagAndProbe/cuts.txt");
  std::string fileSwitches("config/ZeeTagAndProbe/switches.txt");
  
  m_selection = new Selection(fileCuts,fileSwitches);
  configSelection(m_selection);

  // data or MC
  isData_ = m_selection->getSwitch("isData");

  // reading GoodRUN LS 
  std::cout << "[GoodRunLS]::goodRunLS is " << m_selection->getSwitch("goodRunLS") << " isData is " <<  m_selection->getSwitch("isData") << std::endl;

  // to read good run list
  if (m_selection->getSwitch("goodRunLS") && m_selection->getSwitch("isData")) {
    std::string goodRunGiasoneFile = "config/json/goodCollisions2012.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }

  // single electron efficiency
  // simple cuts 2011 (ours and smurfs)
  EgammaCutBasedIDLowPtWPs.push_back("WP95");  // 0
  EgammaCutBasedIDLowPtWPs.push_back("WP90");  // 1
  EgammaCutBasedIDLowPtWPs.push_back("WP85");  // 2
  EgammaCutBasedIDLowPtWPs.push_back("WP80");  // 3
  EgammaCutBasedIDLowPtWPs.push_back("WP70");  // 4
  EgammaCutBasedIDLowPtWPs.push_back("WP70Smurf");  // 5

  EgammaCutBasedIDHighPtWPs.push_back("WP95");  // 0
  EgammaCutBasedIDHighPtWPs.push_back("WP90");  // 1
  EgammaCutBasedIDHighPtWPs.push_back("WP85");  // 2
  EgammaCutBasedIDHighPtWPs.push_back("WP80");  // 3
  EgammaCutBasedIDHighPtWPs.push_back("WP70");  // 4
  EgammaCutBasedIDHighPtWPs.push_back("WP80Smurf");  // 5

  // LH + detector based isolation
  EgammaLHBasedIDLowPtWPs.push_back("LHVeryLooseLowPt");  // 0
  EgammaLHBasedIDLowPtWPs.push_back("LHLooseLowPt");  // 1
  EgammaLHBasedIDLowPtWPs.push_back("LHMediumLowPt");  // 2
  EgammaLHBasedIDLowPtWPs.push_back("LHTightLowPt");  // 3
  EgammaLHBasedIDLowPtWPs.push_back("LHHyperTightLowPt");  // 4

  EgammaLHBasedIDHighPtWPs.push_back("LHVeryLoose");  // 0
  EgammaLHBasedIDHighPtWPs.push_back("LHLoose");  //  1
  EgammaLHBasedIDHighPtWPs.push_back("LHMedium");  // 2
  EgammaLHBasedIDHighPtWPs.push_back("LHTight");  // 3
  EgammaLHBasedIDHighPtWPs.push_back("LHHyperTight");  // 4

  // LH + PFIso
  EgammaLHBasedPFIsoIDLowPtWPs.push_back("LHVeryLoosePFIsoLowPt"); // 0
  EgammaLHBasedPFIsoIDLowPtWPs.push_back("LHLoosePFIsoLowPt");  // 1
  EgammaLHBasedPFIsoIDLowPtWPs.push_back("LHMediumPFIsoLowPt");  // 2
  EgammaLHBasedPFIsoIDLowPtWPs.push_back("LHTightPFIsoLowPt");  // 3
  EgammaLHBasedPFIsoIDLowPtWPs.push_back("LHHyperTightPFIsoLowPt");  // 4

  EgammaLHBasedPFIsoIDHighPtWPs.push_back("LHVeryLoosePFIso");  // 0
  EgammaLHBasedPFIsoIDHighPtWPs.push_back("LHLoosePFIso");  // 1
  EgammaLHBasedPFIsoIDHighPtWPs.push_back("LHMediumPFIso");  // 2
  EgammaLHBasedPFIsoIDHighPtWPs.push_back("LHTightPFIso");  // 3
  EgammaLHBasedPFIsoIDHighPtWPs.push_back("LHHyperTightPFIso");  // 4



  // single electron efficiency
  for (int i=0;i<EgammaCutBasedIDLowPtWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;
      
      char configDir[50];
      sprintf(configDir,"config/%s",EgammaCutBasedIDLowPtWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaCutBasedIDLowPtWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
      //  aSelector.ConfigureEcalCleaner("config/");
      EgammaCutBasedIDLowPt.push_back(aSelector);
    }

  for (int i=0;i<EgammaCutBasedIDHighPtWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;
      
      char configDir[50];
      sprintf(configDir,"config/%s",EgammaCutBasedIDHighPtWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaCutBasedIDHighPtWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
      //  aSelector.ConfigureEcalCleaner("config/");
      EgammaCutBasedIDHighPt.push_back(aSelector);
    }

  for (int i=0;i<EgammaLHBasedIDLowPtWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;

      char configDir[50];
      sprintf(configDir,"config/%s",EgammaLHBasedIDLowPtWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaLHBasedIDLowPtWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
  //  aSelector.ConfigureEcalCleaner("config/");
      EgammaLHBasedIDLowPt.push_back(aSelector);
    }  

  for (int i=0;i<EgammaLHBasedIDHighPtWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;

      char configDir[50];
      sprintf(configDir,"config/%s",EgammaLHBasedIDHighPtWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaLHBasedIDHighPtWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
  //  aSelector.ConfigureEcalCleaner("config/");
      EgammaLHBasedIDHighPt.push_back(aSelector);
    }  

  for (int i=0;i<EgammaLHBasedPFIsoIDLowPtWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;

      char configDir[50];
      sprintf(configDir,"config/%s",EgammaLHBasedPFIsoIDLowPtWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaLHBasedPFIsoIDLowPtWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
  //  aSelector.ConfigureEcalCleaner("config/");
      EgammaLHBasedPFIsoIDLowPt.push_back(aSelector);
    }  

  for (int i=0;i<EgammaLHBasedPFIsoIDHighPtWPs.size();++i)
    {
      CutBasedEleIDSelector aSelector;

      char configDir[50];
      sprintf(configDir,"config/%s",EgammaLHBasedPFIsoIDHighPtWPs[i].c_str());
      std::cout << "===== Configuring " <<  EgammaLHBasedPFIsoIDHighPtWPs[i] << " ElectronID ==========" << std::endl;
      aSelector.ConfigureNoClass(configDir);
  //  aSelector.ConfigureEcalCleaner("config/");
      EgammaLHBasedPFIsoIDHighPt.push_back(aSelector);
    }  

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

  // Reading GoodRUN LS
  std::cout << "[GoodRunLS]::goodRunLS is " << m_selection->getSwitch("goodRunLS") << " isData is " <<  m_selection->getSwitch("isData") << std::endl;
   
  // To read good run list!
  if (m_selection->getSwitch("goodRunLS") && m_selection->getSwitch("isData")) {
    std::string goodRunJsonFile       = "config/json/goodCollisions2012.json";
    setJsonGoodRunList(goodRunJsonFile);
    fillRunLSMap();
  }

}

ZeeTagAndProbe::~ZeeTagAndProbe() { 
  
  delete m_selection;
}


// PDF for probe electrons within the acceptance and loose isolated in the tracker
// the tag is the one with the best match to the Z
// the tag must be within the acceptance, tracker isolated and loose identified
void ZeeTagAndProbe::Loop(const char *treefilesuffix) {
  
  if(fChain == 0) return;
  
  char treename[200];
  sprintf(treename,"%s.root",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesSignal();
  reducedTree.addElectronIdBits();
  reducedTree.addDenominatorFakeBits();
  reducedTree.addIsolations();
  reducedTree.addRunInfos();
  reducedTree.addMore();
  reducedTree.addTrackMomenta();

  // counters
  int allevents   = 0;
  int trigger     = 0;
  int twoele      = 0;
  int eleTot      = 0;
  int eleEta      = 0;
  int elePt       = 0;
  int invmass     = 0;
  int tagId       = 0;
  int tagIsol     = 0;
  int probeIsol   = 0;
  int tagProbeTot = 0;

  // json 
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // triggers
  std::vector<std::string> mask;
  // do one trigger at a time because one need to match the tag with the single Trigger leg
  mask.push_back("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v");
  mask.push_back("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v");
  // triggers for 2012
  mask.push_back("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v");
  mask.push_back("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v");

  // loop over entries
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
    
    // Good Run selection 
    if (isData_ && m_selection->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun  = runNumber;
        lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (isData_ && m_selection->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    allevents++;
    
    // trigger on data
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if (isData_ && !passedHLT ) continue;   
    trigger++;

    // best tag-probe pair = mee closest to Z mass
    float minpull = 100000.;
    float mass    = 1000.;
    float okmass  = 1000.;

    for(int iele1=0; iele1<nEle; iele1++) {
      TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);

      for(int iele2=iele1+1; iele2<nEle; iele2++) {
        TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);
	
	eleTot++;

        if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[iele1]) || !m_selection->passCut("etaEleAcc",etaEle[iele2]) ) ) continue;
	eleEta++;

        if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc",electron1.Pt()) || !m_selection->passCut("ptEleAcc",electron2.Pt()) ) ) continue;
	elePt++;

        
	
        mass = (electron1+electron2).M();
        float pull=fabs(mass-91.1876);
	if(pull < minpull) {
	  okmass  = mass;
          minpull = pull;
          electrons[0] = iele1;
          electrons[1] = iele2;
        }
      }
    }

    if (okmass<999) twoele++;

    // start the tag & probe
    if( m_selection->passCut("meeWindow",okmass) ) {

      if ( okmass>120 || okmass<60 ) cout << "BACO!" << endl;
      invmass++;
      
      for(int iele=0; iele<2; ++iele) {
        int tag=electrons[0], probe=electrons[1];
        if(iele==1) {
          tag=electrons[1]; 
          probe=electrons[0];
        }

        bool mcmatch=false;
        if(!isData_) mcmatch = mcMatches(probe,tag);

        // match the tag electron with HLT candidate
        bool tagMatch = triggerMatch(etaEle[tag],phiEle[tag],0.2);
        bool tagProbe = triggerMatch(etaEle[probe],phiEle[probe],0.2);
        //        if((!tagMatch) || (!tagProbe)) continue;

        TLorentzVector probeP4(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);
        TLorentzVector tagP4(pxEle[tag],    pyEle[tag],  pzEle[tag],  energyEle[tag]);

	// various about probe
        int charge = chargeEle[probe];
        float pt   = probeP4.Pt();
        float eta  = etaEle[probe];

        // apply the electron ID tight (WP70Smurf) on the tag electron
        bool tagIdentified, tagIsolated, tagConvRej;
        tagIdentified = tagIsolated = tagConvRej = false;
        isEleID(&EgammaCutBasedIDLowPt[5],tag,&tagIdentified,&tagIsolated,&tagConvRej);
        // override: patch
        tagIdentified = tagIsolated = tagConvRej = true;
        if (tagIdentified && tagIsolated) {

          // look at the probe
          int CutBasedId[6], CutBasedIdOnlyID[6], CutBasedIdOnlyIso[6], CutBasedIdOnlyConv[6];
          int LHBasedId[5], LHBasedIdOnlyID[5], LHBasedIdOnlyIso[5], LHBasedIdOnlyConv[5];
          int LHBasedPFIsoId[5], LHBasedPFIsoIdOnlyID[5], LHBasedPFIsoIdOnlyIso[5], LHBasedPFIsoIdOnlyConv[5];
          
          int cutssize = (pt > 20) ? EgammaCutBasedIDHighPtWPs.size() : EgammaCutBasedIDLowPtWPs.size();

          for (int icut=0;icut<cutssize;++icut)
          {
            bool isEleIDCutBased, isIsolCutBased, isConvRejCutBased;
            isEleIDCutBased = isIsolCutBased = isConvRejCutBased = false;
            if(pt > 20) isEleID(&EgammaCutBasedIDHighPt[icut],probe,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
            else {
              isEleID(&EgammaCutBasedIDLowPt[icut],probe,&isEleIDCutBased,&isIsolCutBased,&isConvRejCutBased);
              if(TString(EgammaCutBasedIDLowPtWPs[icut].c_str()).Contains("Smurf") && pt<15.) {
                isEleIDCutBased = isEleIDCutBased && (fbremEle[probe]>0.15 || (fabs(etaEle[probe])<1.0 && eSuperClusterOverPEle[probe]>0.95));
              }
            }
            CutBasedIdOnlyID[icut] = (isEleIDCutBased) ? 1 : 0;
            CutBasedIdOnlyIso[icut] = (isIsolCutBased) ? 1 : 0;
            CutBasedIdOnlyConv[icut] = (isConvRejCutBased) ? 1 : 0;
            CutBasedId[icut] = (isEleIDCutBased && isIsolCutBased && isConvRejCutBased) ? 1 : 0;
          }
          
          cutssize = (pt > 20) ? EgammaLHBasedIDHighPtWPs.size() : EgammaLHBasedIDLowPtWPs.size();
          
          for (int icut=0;icut<cutssize;++icut)
            {
              bool isEleIDLHBased, isIsolLHBased, isConvRejLHBased;
              isEleIDLHBased = isIsolLHBased = isConvRejLHBased = false;
              if(pt > 20) isEleID(&EgammaLHBasedIDHighPt[icut],probe,&isEleIDLHBased,&isIsolLHBased,&isConvRejLHBased);
              else isEleID(&EgammaLHBasedIDLowPt[icut],probe,&isEleIDLHBased,&isIsolLHBased,&isConvRejLHBased);
              LHBasedIdOnlyID[icut] = (isEleIDLHBased) ? 1 : 0;
              LHBasedIdOnlyIso[icut] = (isIsolLHBased) ? 1 : 0;
              LHBasedIdOnlyConv[icut] = (isConvRejLHBased) ? 1 : 0;
              LHBasedId[icut] = (isEleIDLHBased && isIsolLHBased && isConvRejLHBased) ? 1 : 0;
            }

          cutssize = (pt > 20) ? EgammaLHBasedPFIsoIDHighPtWPs.size() : EgammaLHBasedPFIsoIDLowPtWPs.size();
          
          for (int icut=0;icut<cutssize;++icut)
            {
              bool isEleIDLHBasedPFIso, isIsolLHBasedPFIso, isConvRejLHBasedPFIso;
              isEleIDLHBasedPFIso = isIsolLHBasedPFIso = isConvRejLHBasedPFIso = false;
              if(pt > 20) isEleID(&EgammaLHBasedPFIsoIDHighPt[icut],probe,&isEleIDLHBasedPFIso,&isIsolLHBasedPFIso,&isConvRejLHBasedPFIso);
              else isEleID(&EgammaLHBasedPFIsoIDLowPt[icut],probe,&isEleIDLHBasedPFIso,&isIsolLHBasedPFIso,&isConvRejLHBasedPFIso);
              LHBasedPFIsoIdOnlyID[icut] = (isEleIDLHBasedPFIso) ? 1 : 0;
              LHBasedPFIsoIdOnlyIso[icut] = (isIsolLHBasedPFIso) ? 1 : 0;
              LHBasedPFIsoIdOnlyConv[icut] = (isConvRejLHBasedPFIso) ? 1 : 0;
              LHBasedPFIsoId[icut] = (isEleIDLHBasedPFIso && isIsolLHBasedPFIso && isConvRejLHBasedPFIso) ? 1 : 0;
            }
        
          // some eleID variables
          float HoE, s1s9, s9s25, phiwidth, etawidth, deta, dphi, fbrem, see, spp, eleopout, eopout, eop, eseedopin, nbrems, recoFlag, EleSCEta, EleSCPhi;
          float EleSCEt,EtaSeed,PhiSeed,ESeed,IEtaSeed,IPhiSeed,EtaCrySeed,PhiCrySeed,IEtaCrySeed,IPhiCrySeed;
          float oneoveremoneoverp, eledeta, d0, ip3d, ip3ds, kfhits, kflayers, kfchi2, e1x5e5x5, dcot, dist;
          float detacalo, dphicalo, sep, dz, gsfchi2, emax, etop, ebottom, eleft, eright,
            e2nd, e2x5right, e2x5left, e2x5top, e2x5bottom, 
            e2x5max, e1x5, e2x2, e3x3, e5x5, r9, nclu,
            phi, scenergy, scrawenergy, scesenergy;

          phi  = phiEle[probe];

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
              if(seedclu!=-1) { // may happen because PFBC from photons are saved...
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
                EtaSeed=-999;
                PhiSeed=-999;
                ESeed=-999.;
                IEtaSeed=-999.;
                IPhiSeed=-999.;
                EtaCrySeed=-999.;
                PhiCrySeed=-999.;
                IEtaCrySeed=-999.;
                IPhiCrySeed=-999.;
              }
            } else {
              s9s25 = 999.;
              see = 999.;
              spp = 999.;
            }
          }

          // this only happens in 42X EE because of a bug in the EcalClusterTools
          if(isnan(see) || isnan(spp) || isnan(sep)) continue;

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
          hwwbdts[1] = eleBDTWithIso(fMVAHWWWithIso,probe);
          float hzzbdts[4];
          hzzbdts[0] = eleBDT(fMVAHZZDanV0,probe);
          hzzbdts[1] = eleBDT(fMVAHZZSiV0,probe);
          hzzbdts[2] = eleBDT(fMVAHZZSiV1,probe);
          //hzzbdts[3] = eleBDT(fMVAHZZSiDanV2,probe);
          hzzbdts[3] = mvaidnontrigEle[probe];
          float newhwwbdts[4];
          newhwwbdts[0] = eleBDT(fMVAHWWDanV0,probe);
          newhwwbdts[1] = eleBDT(fMVAHWWSiV0,probe);
          newhwwbdts[2] = eleBDT(fMVAHWWSiV1,probe);
          //newhwwbdts[3] = eleBDT(fMVAHWWSiDanV2,probe);
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
          
          float GenEnergy = genEleEnergy(probe);
          float GenEnergyStatus1 = genEleEnergy(probe,1);
          float GenEnergyStatus3 = genEleEnergy(probe,3);

          // fill the reduced tree
	  reducedTree.fillVariables(eleopout,eopout,eop,HoE,deta,dphi,s9s25,s1s9,see,spp,fbrem,
                                    nbrems,misshits,dcot,dist,pt,eta,charge,phiwidth,etawidth,
                                    oneoveremoneoverp,eledeta,d0,ip3d,ip3ds,kfhits,kflayers,kfchi2,e1x5e5x5,ecalseed,matchConv,
                                    isEleEB,isEleEE);
          reducedTree.fillVariables2(detacalo, dphicalo, sep, dz, gsfchi2, emax, etop, ebottom, eleft, eright,
                                     e2nd, e2x5right, e2x5left, e2x5top, e2x5bottom, 
                                     e2x5max, e1x5, e2x2, e3x3, e5x5, r9, nclu,
                                     phi, scenergy, scrawenergy, scesenergy,eseedopin);
          reducedTree.fillCluterInfos(EleSCEt,EleSCEta,EleSCPhi,EtaSeed,PhiSeed,ESeed,IEtaSeed,IPhiSeed,EtaCrySeed,PhiCrySeed,IEtaCrySeed,IPhiCrySeed);
          reducedTree.fillAttributesSignal(okmass,1,GenEnergy,GenEnergyStatus1,GenEnergyStatus3);
          reducedTree.fillIsolations(trkIso,ecalIso,hcalIso,
                                     pfCombinedIsoEle[probe],
                                     chaPfIso, neuPfIso, phoPfIso);
          reducedTree.fillMore(nPV,rhoFastjet,hwwbdts,newhwwbdts,hzzbdts,pfmva,lh);
	  reducedTree.fillTrackMomenta(pcomb,pmodegsf,pmeangsf,pmeankf,pterrorgsf,pterrorkf,perrorele);
          reducedTree.fillCutBasedIDBits(CutBasedId,CutBasedIdOnlyID,CutBasedIdOnlyIso,CutBasedIdOnlyConv);
          reducedTree.fillLHBasedIDBits(LHBasedId,LHBasedIdOnlyID,LHBasedIdOnlyIso,LHBasedIdOnlyConv);
          reducedTree.fillLHBasedPFIsoIDBits(LHBasedPFIsoId,LHBasedPFIsoIdOnlyID,LHBasedPFIsoIdOnlyIso,LHBasedPFIsoIdOnlyConv);
          reducedTree.fillFakeRateDenomBits(-1.,isDenomFake(probe),isDenomFake_smurfs(probe));
          reducedTree.fillBDTBasedIDBits(passEleBDT(pt,EleSCEta,hwwbdts[0]));
	  reducedTree.fillCiCBasedIDBits(cic);
          reducedTree.fillRunInfos(runNumber, lumiBlock, eventNumber, nPU, mcmatch);
          reducedTree.store();
        } // tag identified

      } // loop over the 2 Z electrons
      
    } // end tag and probe
    
  } // loop over events
  
  cout << "statistics from Tag and Probe: " << endl;
  cout << "allevents   = " << allevents << endl;
  cout << "trigger     = " << trigger << endl;
  cout << "twoele      = " << twoele  << endl;
  cout << "invmass     = " << invmass << endl;
  cout << "tagProbeTot = " << tagProbeTot << endl;
  cout << "tagId       = " << tagId       << endl;
  cout << "tagIsol     = " << tagIsol     << endl;
  cout << "probeIsol   = " << probeIsol   << endl;  
  cout << "statistics from Tag and Probe - electrons: " << endl;
  cout << "eleTot      = " << eleTot  << endl;
  cout << "eleEta      = " << eleEta  << endl;
  cout << "elePt       = " << elePt   << endl;

  reducedTree.save();
}


void ZeeTagAndProbe::isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

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
  //  selector->SetEgammaCutBasedID ( anaUtils.electronIdVal(eleIdCutsEle[eleIndex],eleIdLoose) );
  selector->SetLikelihood( lh );
  selector->SetEcalIsolation( (dr03EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  selector->SetTrkIsolation( (dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  selector->SetHcalIsolation( (dr03HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3)/pEle.Pt() );
  float combinedIso = 0.0;
  if (isEleEB) combinedIso = dr03TkSumPtEle[eleIndex] + TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + dr03HcalTowerSumEtFullConeEle[eleIndex];
  else combinedIso = dr03TkSumPtEle[eleIndex] + dr03EcalRecHitSumEtEle[eleIndex] + dr03HcalTowerSumEtFullConeEle[eleIndex];
  selector->SetCombinedIsolation( (combinedIso - rhoFastjet*TMath::Pi()*0.3*0.3) / pEle.Pt() );

  selector->SetCombinedPFIsolation( (pfCombinedIsoEle[eleIndex]) / pEle.Pt() );
  // selector->SetCombinedPFIsolation( (pfCombinedIsoEle[eleIndex] - rhoFastjet*TMath::Pi()*0.4*0.4) / pEle.Pt() );

  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );
  selector->SetHasMatchedConversion ( hasMatchedConversionEle[eleIndex] );

  // ECAL cleaning variables
  //selector->m_cleaner->SetE1(e1);
  //selector->m_cleaner->SetE4SwissCross(e4SwissCross);
  //selector->m_cleaner->SetFiducialFlag(fidFlagSC);
  //selector->m_cleaner->SetSeedFlag(seedRecHitFlag);
  //selector->m_cleaner->SetSeedTime(seedTime);
  //selector->m_cleaner->SetSeedChi2(seedChi2);

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputNoClassEleId();
  *isolOutput = selector->outputNoClassIso();
  *convRejOutput = selector->outputNoClassConv();

}

void ZeeTagAndProbe::isEleID(CiCBasedEleSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];

  TVector3 pEle(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eseedopin, eop, eseed,  sce, scet,sceta, rawe;
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
    eseed = seedEnergySC[sc];
    sce = energySC[sc];
    rawe = rawEnergySC[sc];
    scet = energySC[sc]/TMath::CosH(etaSC[sc]);
    sceta = etaSC[sc];
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
      eseed = seedEnergyPFSC[sc];
      sce = energyPFSC[sc];
      rawe = rawEnergyPFSC[sc];
      scet = energyPFSC[sc]/TMath::CosH(etaPFSC[sc]);
      sceta = etaPFSC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  
  eseedopin = eSeedOverPoutEle[eleIndex]*(1-fbremEle[eleIndex]);
  //  std::cout << eop << "," <<  eSeedOverPoutEle[eleIndex] << ","<< eseedopin*sce/rawe << "," << fbremEle[eleIndex] << "," << nbremsEle[eleIndex] << std::endl;;
  
  //  float lh=likelihoodRatio(eleIndex,*LH);
  float lh = eleIdLikelihoodEle[eleIndex];
  selector->reset();
  selector->SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  selector->SetRecoFlag(recoFlagsEle[eleIndex]);
  selector->SetSCEt(scet);
  selector->SetSCEta(sceta);
  selector->SetHOverE( HoE );
  selector->SetDEta( deta );
  selector->SetDPhiIn( dphiin );
  selector->SetBremFraction( fbrem );
  selector->SetSigmaEtaEta( see );
  selector->SetEOverPin( eop );
  selector->SetESeedOverPin( eseedopin );
  selector->SetEcalIsolation( dr04EcalRecHitSumEtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.4*0.4 );
  selector->SetTrkIsolation( dr03TkSumPtEle[eleIndex] - rhoFastjet*TMath::Pi()*0.3*0.3 );
  selector->SetHcalIsolation( dr04HcalTowerSumEtFullConeEle[eleIndex] - rhoFastjet*TMath::Pi()*0.4*0.4);
  selector->SetMissingHits( expInnerLayersGsfTrack[gsf] );
  selector->SetConvDist( fabs(convDistEle[eleIndex]) );
  selector->SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
//   selector->m_cleaner->SetE1(e1);
//   selector->m_cleaner->SetE4SwissCross(e4SwissCross);
//   selector->m_cleaner->SetFiducialFlag(fidFlagSC);
//   selector->m_cleaner->SetSeedFlag(seedRecHitFlag);
//   selector->m_cleaner->SetSeedTime(seedTime);
//   selector->m_cleaner->SetSeedChi2(seedChi2);

  //  return selector->output(); // class dependent result
  *eleIdOutput = selector->outputEleId();
  *isolOutput = selector->outputIso();
  *convRejOutput = selector->outputConv();

}

// denominator for fake rate: for HtoWW, egamma triggers
int ZeeTagAndProbe::isDenomFake(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);

  // match with the HLT firing candidates
  //  bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
  //  if (!HLTmatch) isGoodDenom = false;
  
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
  return 0;
}

// denominator for fake rate: for HtoWW, egamma triggers, same as smurfs
int ZeeTagAndProbe::isDenomFake_smurfs(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // match with the HLT firing candidates
  // bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
  // if (!HLTmatch) isGoodDenom = false;
  
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
  float dzEle  = PVzPV[0] - trackVzGsfTrack[gsf];
  if (fabs(dxyEle)>0.02) isGoodDenom = false;
  if (fabs(dzEle)>0.10)  isGoodDenom = false;

  if(isGoodDenom) return 1;
  return 0;
}


void ZeeTagAndProbe::configSelection(Selection* selection) {

  m_selection->addSwitch("isData");
  m_selection->addSwitch("goodRunLS");
  m_selection->addCut("meeWindow");
  m_selection->addCut("etaEleAcc");
  m_selection->addCut("ptEleAcc");
  m_selection->summary();

}


bool ZeeTagAndProbe::mcMatches(int probe, int tag) {
  
  bool probematch=true;
  bool tagmatch=true;

  int ep,em;
  ep=em=-1;
  
  for(int imc=0; imc<16; ++imc) {
    if(idMc[imc]==11 && abs(idMc[mothMc[imc]])==23) em=imc;
    else if(idMc[imc]==-11 && abs(idMc[mothMc[imc]])==23) ep=imc;
    if(em>-1 && ep>-1) break;
  }
  
  TVector3 probeP(pxEle[probe],pyEle[probe],pzEle[probe]);
  int probemc, tagmc;
  probemc=tagmc=-1;
  if(chargeEle[probe]>0) {
    probemc=ep;
    tagmc=em;
  } else {
    probemc=em;
    tagmc=ep;
  }
  TVector3 probemcP;
  probemcP.SetMagThetaPhi(pMc[probemc],thetaMc[probemc],phiMc[probemc]);
  if(probemcP.DeltaR(probeP)>0.1) probematch=false;
  
  if(tag>-1) {
    TVector3 tagP(pxEle[tag],pyEle[tag],pzEle[tag]);
    TVector3 tagmcP;
    tagmcP.SetMagThetaPhi(pMc[tagmc],thetaMc[tagmc],phiMc[tagmc]);
    if(tagmcP.DeltaR(tagP)>0.1) tagmatch=false;
  }
  return (probematch && tagmatch);
}

float ZeeTagAndProbe::genEleEnergy(int ele, int status, float deltar) {
  int gen1,gen2;
  gen1=gen2=-1;
  for(int imc=0; imc<20; ++imc) {
    // mother or grandmother should be the Z (to allow the match with status 1 or 3 electron
    if(abs(idMc[imc])==11 && (abs(idMc[mothMc[imc]])==23 || abs(idMc[mothMc[mothMc[imc]]]))) {
      if(status>0 && statusMc[imc]!=status) continue;
      if(gen1==-1) gen1=imc;
      else gen2=imc;
    }
  }
  if(gen1<0 || gen2<0) return -1.0;

  TVector3 eleP(pxEle[ele],pyEle[ele],pzEle[ele]);
  TVector3 mcele1P, mcele2P;
  mcele1P.SetMagThetaPhi(pMc[gen1],thetaMc[gen1],phiMc[gen1]);
  mcele2P.SetMagThetaPhi(pMc[gen2],thetaMc[gen2],phiMc[gen2]);
  float dr1=mcele1P.DeltaR(eleP);
  float dr2=mcele2P.DeltaR(eleP);
  float mindr=TMath::Min(dr1,dr2);
  if(mindr>deltar) return -2.0; 

  if(dr1<dr2) return mcele1P.Mag();
  else return mcele2P.Mag();
}
