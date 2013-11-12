#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/Utils.hh"
#include <iostream>
#include <math.h>
#include <string.h>

using namespace bits;

CutBasedEleIDSelector::CutBasedEleIDSelector() {

  // if the variable value is not initialized, it is not used
  m_useHOverE = false;
  m_useS9S25 = false;
  m_useDEta = false;
  m_useDPhiIn = false;
  m_useDPhiOut = false;
  m_useInvEminusInvP = false;
  m_useBremFraction = false;
  m_useSigmaEtaEta = false;
  m_useSigmaPhiPhi = false;
  m_useEOverPout = false;
  m_useEOverPin = false;
  m_useEcalIso = false;
  m_useTrkIso = false;
  m_useHcalIso = false;
  m_useCombIso = false;
  m_useCombPFIso = false;
  m_useDistConv = false;
  m_useDcotConv = false;
  m_useMatchedConv = false;
  m_useMissingHits = false;
  m_useLikelihood = false;
  m_electronClassInitialised = false;
  m_egammaCutBasedInitialised = false;

  m_applyIDOnPFlow = true;

  m_doEcalCleaning = false;

  // fiducial flag HAS to be set, otherwise the selector doesn't know which slection apply (EB or EE)
  m_fiducialflag = -1;
  
  // if true, the selector applies the same ID for PFlow and standard electrons
  m_applyIDOnPFlow = true;
  // if not set, assign the flag to ECAL driven
  m_recoflag = 2;

}

CutBasedEleIDSelector::~CutBasedEleIDSelector() {}

void CutBasedEleIDSelector::Configure(const char *configDir) {

  m_classDep = true;

  std::string fileSwitches = std::string(configDir) + "/electronIDSwitches.txt";

  std::string goldenCutsEB = std::string(configDir) + "/electronIDGoldenCutsEB.txt";
  std::string bigbremCutsEB = std::string(configDir) + "/electronIDBigBremCutsEB.txt";
  std::string narrowCutsEB = std::string(configDir) + "/electronIDNarrowCutsEB.txt";
  std::string showeringCutsEB = std::string(configDir) + "/electronIDShoweringCutsEB.txt";

  std::string goldenCutsEE = std::string(configDir) + "/electronIDGoldenCutsEE.txt";
  std::string bigbremCutsEE = std::string(configDir) + "/electronIDBigBremCutsEE.txt";
  std::string narrowCutsEE = std::string(configDir) + "/electronIDNarrowCutsEE.txt";
  std::string showeringCutsEE = std::string(configDir) + "/electronIDShoweringCutsEE.txt";
  
  Selection *goldenSelectionEB = new Selection(goldenCutsEB, fileSwitches);
  Selection *bigbremSelectionEB = new Selection(bigbremCutsEB, fileSwitches);
  Selection *narrowSelectionEB = new Selection(narrowCutsEB, fileSwitches);
  Selection *showeringSelectionEB = new Selection(showeringCutsEB, fileSwitches);

  Selection *goldenSelectionEE = new Selection(goldenCutsEE, fileSwitches);
  Selection *bigbremSelectionEE = new Selection(bigbremCutsEE, fileSwitches);
  Selection *narrowSelectionEE = new Selection(narrowCutsEE, fileSwitches);
  Selection *showeringSelectionEE = new Selection(showeringCutsEE, fileSwitches);

  m_EgammaCutIDSelection.push_back(goldenSelectionEB);
  m_EgammaCutIDSelection.push_back(bigbremSelectionEB);
  m_EgammaCutIDSelection.push_back(narrowSelectionEB);
  m_EgammaCutIDSelection.push_back(showeringSelectionEB);
  m_EgammaCutIDSelection.push_back(goldenSelectionEE);
  m_EgammaCutIDSelection.push_back(bigbremSelectionEE);
  m_EgammaCutIDSelection.push_back(narrowSelectionEE);
  m_EgammaCutIDSelection.push_back(showeringSelectionEE);

  std::vector<Selection*>::iterator selectionItr;
  for(selectionItr=m_EgammaCutIDSelection.begin();selectionItr!=m_EgammaCutIDSelection.end();selectionItr++) {
    Selection *eleSelection = *selectionItr;
    eleSelection->addSwitch("egammaCutBased");
    eleSelection->addCut("hOverE");
    eleSelection->addCut("s9s25");
    eleSelection->addCut("deta");
    eleSelection->addCut("dphiIn");
    eleSelection->addCut("dphiOut");
    eleSelection->addCut("invEminusInvP");
    eleSelection->addCut("bremFraction"); 
    eleSelection->addCut("sigmaEtaEta");
    eleSelection->addCut("sigmaPhiPhi");
    eleSelection->addCut("eOverPout");
    eleSelection->addCut("eOverPin");
    eleSelection->addCut("ecalIso");
    eleSelection->addCut("trackerIso");
    eleSelection->addCut("hcalIso");
    eleSelection->addCut("combIso");
    eleSelection->addCut("combPFIso");
    eleSelection->addCut("missHits");
    eleSelection->addCut("distConv");
    eleSelection->addCut("dcotConv");
    eleSelection->addSwitch("matchedConv");
    eleSelection->addCut("likelihood");
    eleSelection->summary();
  }


  m_electronCounter.SetTitle("SINGLE ELECTRON COUNTER");
  m_electronCounter.AddVar("electrons");
  m_electronCounter.AddVar("electronsOnlyID");
  m_electronCounter.AddVar("electronsOnlyIso");
  m_electronCounter.AddVar("electronsOnlyConv");
  m_electronCounter.AddVar("egammaCutBased");
  m_electronCounter.AddVar("hOverE");
  m_electronCounter.AddVar("s9s25");
  m_electronCounter.AddVar("deta");
  m_electronCounter.AddVar("dphiIn");
  m_electronCounter.AddVar("dphiOut");
  m_electronCounter.AddVar("invEminusInvP");
  m_electronCounter.AddVar("bremFraction");
  m_electronCounter.AddVar("sigmaEtaEta");
  m_electronCounter.AddVar("sigmaPhiPhi");
  m_electronCounter.AddVar("eOverPout");
  m_electronCounter.AddVar("eOverPin");
  m_electronCounter.AddVar("ecalCleaning");
  m_electronCounter.AddVar("trackerIso");
  m_electronCounter.AddVar("hcalIso");
  m_electronCounter.AddVar("combIso");
  m_electronCounter.AddVar("combPFIso");
  m_electronCounter.AddVar("missHits");
  m_electronCounter.AddVar("vetoConv");
  m_electronCounter.AddVar("matchedConv");
  m_electronCounter.AddVar("finalCustomEleIDOnlyID");
  m_electronCounter.AddVar("finalCustomEleIDOnlyIso");
  m_electronCounter.AddVar("finalCustomEleIDOnlyConv");
  m_electronCounter.AddVar("likelihood");


}

void CutBasedEleIDSelector::ConfigureNoClass(const char *configDir) 
{

  m_classDep = false;

  std::string fileSwitches = std::string(configDir) + "/electronIDSwitches.txt";

  std::string cutsEB = std::string(configDir) + "/electronIDCutsEB.txt";
  std::string cutsEE = std::string(configDir) + "/electronIDCutsEE.txt";
  
  Selection *selectionEB = new Selection(cutsEB, fileSwitches);
  Selection *selectionEE = new Selection(cutsEE, fileSwitches);

  m_EgammaCutIDSelection.push_back(selectionEB);
  m_EgammaCutIDSelection.push_back(selectionEE);

  std::vector<Selection*>::iterator selectionItr;
  for(selectionItr=m_EgammaCutIDSelection.begin();selectionItr!=m_EgammaCutIDSelection.end();selectionItr++) {
    Selection *eleSelection = *selectionItr;
    eleSelection->addSwitch("egammaCutBased");
    eleSelection->addCut("hOverE");
    eleSelection->addCut("s9s25");
    eleSelection->addCut("deta");
    eleSelection->addCut("dphiIn");
    eleSelection->addCut("dphiOut");
    eleSelection->addCut("invEminusInvP");
    eleSelection->addCut("bremFraction"); 
    eleSelection->addCut("sigmaEtaEta");
    eleSelection->addCut("sigmaPhiPhi");
    eleSelection->addCut("eOverPout");
    eleSelection->addCut("eOverPin");
    eleSelection->addCut("ecalIso");
    eleSelection->addCut("trackerIso");
    eleSelection->addCut("hcalIso");
    eleSelection->addCut("combIso");
    eleSelection->addCut("combPFIso");
    eleSelection->addCut("missHits");
    eleSelection->addCut("distConv");
    eleSelection->addCut("dcotConv");
    eleSelection->addSwitch("matchedConv");
    eleSelection->addCut("likelihood0");
    eleSelection->addCut("likelihood1");
    eleSelection->summary();
  }


  m_electronCounter.SetTitle("SINGLE ELECTRON COUNTER");
  m_electronCounter.AddVar("electrons");
  m_electronCounter.AddVar("electronsOnlyID");
  m_electronCounter.AddVar("electronsOnlyIso");
  m_electronCounter.AddVar("electronsOnlyConv");
  m_electronCounter.AddVar("egammaCutBased");
  m_electronCounter.AddVar("hOverE");
  m_electronCounter.AddVar("s9s25");
  m_electronCounter.AddVar("deta");
  m_electronCounter.AddVar("dphiIn");
  m_electronCounter.AddVar("dphiOut");
  m_electronCounter.AddVar("invEminusInvP");
  m_electronCounter.AddVar("bremFraction");
  m_electronCounter.AddVar("sigmaEtaEta");
  m_electronCounter.AddVar("sigmaPhiPhi");
  m_electronCounter.AddVar("eOverPout");
  m_electronCounter.AddVar("eOverPin");
  m_electronCounter.AddVar("ecalCleaning");
  m_electronCounter.AddVar("ecalIso");
  m_electronCounter.AddVar("trackerIso");
  m_electronCounter.AddVar("hcalIso");
  m_electronCounter.AddVar("combIso");
  m_electronCounter.AddVar("combPFIso");
  m_electronCounter.AddVar("missHits");
  m_electronCounter.AddVar("vetoConv");
  m_electronCounter.AddVar("matchedConv");
  m_electronCounter.AddVar("likelihood");
  m_electronCounter.AddVar("finalCustomEleID");
  m_electronCounter.AddVar("finalCustomEleIDOnlyID");
  m_electronCounter.AddVar("finalCustomEleIDOnlyIso");
  m_electronCounter.AddVar("finalCustomEleIDOnlyConv");


}

void CutBasedEleIDSelector::ConfigureEcalCleaner(const char *configDir) {

  m_doEcalCleaning = true;

  m_cleaner = new EcalCleaner();
  m_cleaner->Configure(configDir);

}

Selection* CutBasedEleIDSelector::GetSelectionNoClass(const char *EcalSubDet) {

  if(m_EgammaCutIDSelection.size() != 2) {
    std::cout << "CutBasedEleIDSelector::GetSelectionNoClass ERROR! The selection seems to be configured with classes" << std::endl;
    return 0;
  }
  if(strcmp(EcalSubDet,"EB") == 0) return m_EgammaCutIDSelection[0];
  else if(strcmp(EcalSubDet,"EE") == 0) return m_EgammaCutIDSelection[1];

  std::cout << "CutBasedEleIDSelector::GetSelectionNoClass ERROR! Provide subdet string = EB / EE" << std::endl; 
  return 0;

}

bool CutBasedEleIDSelector::output() {
  
  if ( m_fiducialflag == -1 ) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CutBasedEleIDSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return false;
  }
  
  if ( !m_classDep ) {
    cout << "CutBasedEleIDSelector::output() method is intended for class based ele ID, while you configured for class independent ID. Use outputNoClass() instead." << endl;
    return false;
  }

  m_electronCounter.IncrVar("electrons");

  Utils anaUtils;
  bool isPFlowOnly = !anaUtils.electronRecoType(m_recoflag, isEcalDriven);
  if ( isPFlowOnly && !m_applyIDOnPFlow ) return true;

  if (!outputEleId())
    return false;
  if (!outputIso())
    return false;
  if (!outputConv())
    return false;
  m_electronCounter.IncrVar("finalCustomEleID");

  return true;  

}
bool CutBasedEleIDSelector::outputEleId() {

  if ( m_fiducialflag == -1 ) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CutBasedEleIDSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return false;
  }

  if ( !m_classDep ) {
    cout << "CutBasedEleIDSelector::output() method is intended for class based ele ID, while you configured for class independent ID. Use outputNoClass() instead." << endl;
    return false;
  }


  m_electronCounter.IncrVar("electronsOnlyID");

  Utils anaUtils;
  bool isPFlowOnly = !anaUtils.electronRecoType(m_recoflag, isEcalDriven);
  if ( isPFlowOnly && !m_applyIDOnPFlow ) return true;

  bool eleInGap = anaUtils.fiducialFlagECAL(m_fiducialflag, isEBEEGap);

  if ( eleInGap ) return false;

  int GsfClass = m_electronClass;
  int offset=0;
  bool eeElectron = anaUtils.fiducialFlagECAL(m_fiducialflag, isEE);
  if(eeElectron) {
    offset=4;
  }

  Selection *selection;
  if(GsfClass==GOLDEN) 
    selection=m_EgammaCutIDSelection[offset];
  else if(GsfClass==BIGBREM)
    selection=m_EgammaCutIDSelection[offset+1];
  else if(GsfClass==NARROW)
    selection=m_EgammaCutIDSelection[offset+2];
  else if(GsfClass==SHOWERING) {
    selection=m_EgammaCutIDSelection[offset+3];
  }


  if(selection->getSwitch("egammaCutBased") && 
     m_egammaCutBased ) m_electronCounter.IncrVar("egammaCutBased");

  if(selection->getSwitch("likelihood") && 
     !selection->passCut("likelihood", fabs(m_Likelihood) )) return false;
  m_electronCounter.IncrVar("likelihood");

  if(selection->getSwitch("hOverE") && 
     !selection->passCut("hOverE", fabs(m_HOverE) )) return false;  
  m_electronCounter.IncrVar("hOverE");
  
  if(selection->getSwitch("s9s25") && 
     !selection->passCut("s9s25", m_S9S25)) return false; 
  m_electronCounter.IncrVar("s9s25");
  
  if(selection->getSwitch("deta") && 
     !selection->passCut("deta", fabs(m_DEta) )) return false; 
  m_electronCounter.IncrVar("deta");

  if(selection->getSwitch("dphiIn") && 
     !selection->passCut("dphiIn", fabs(m_DPhiIn) )) return false; 
  m_electronCounter.IncrVar("dphiIn");

  if(selection->getSwitch("dphiOut") && 
     !selection->passCut("dphiOut", fabs(m_DPhiOut) )) return false; 
  m_electronCounter.IncrVar("dphiOut");
    
  if(selection->getSwitch("invEminusInvP") && 
     !selection->passCut("invEminusInvP", m_InvEminusInvP)) return false; 
  m_electronCounter.IncrVar("invEminusInvP");
    
  if(selection->getSwitch("bremFraction") && 
     !selection->passCut("bremFraction", m_BremFraction)) return false; 
  m_electronCounter.IncrVar("bremFraction");
  
  if(selection->getSwitch("sigmaEtaEta") && 
     !selection->passCut("sigmaEtaEta", m_SigmaEtaEta)) return false; 
  m_electronCounter.IncrVar("sigmaEtaEta");

  if(selection->getSwitch("sigmaPhiPhi") && 
     !selection->passCut("sigmaPhiPhi", m_SigmaPhiPhi)) return false; 
  m_electronCounter.IncrVar("sigmaPhiPhi");

  if(selection->getSwitch("eOverPout") && 
     !selection->passCut("eOverPout", m_EOverPout)) return false; 
  m_electronCounter.IncrVar("eOverPout");

  if(selection->getSwitch("eOverPin") && 
     !selection->passCut("eOverPin", m_EOverPin)) return false; 
  m_electronCounter.IncrVar("eOverPin");

  if(m_doEcalCleaning && !m_cleaner->output()) return false;
  m_electronCounter.IncrVar("ecalCleaning");
  m_electronCounter.IncrVar("finalCustomEleIDOnlyID");

  return true;


}


bool CutBasedEleIDSelector::outputIso() 
{
  if ( m_fiducialflag == -1 ) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CutBasedEleIDSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return false;
  }

  if ( !m_classDep ) {
    cout << "CutBasedEleIDSelector::output() method is intended for class based ele ID, while you configured for class independent ID. Use outputNoClass() instead." << endl;
    return false;
  }
  
    int GsfClass = m_electronClass;
  int offset=0;
  Utils anaUtils;
  bool eeElectron = anaUtils.fiducialFlagECAL(m_fiducialflag, isEE);
  if(eeElectron) {
    offset=4;
  }

  Selection *selection;
  if(GsfClass==GOLDEN) 
    selection=m_EgammaCutIDSelection[offset];
  else if(GsfClass==BIGBREM)
    selection=m_EgammaCutIDSelection[offset+1];
  else if(GsfClass==NARROW)
    selection=m_EgammaCutIDSelection[offset+2];
  else if(GsfClass==SHOWERING) {
    selection=m_EgammaCutIDSelection[offset+3];
  }

  m_electronCounter.IncrVar("electronsOnlyIso");

  if(!selection->getSwitch("combIso") && !selection->getSwitch("combPFIso")) {
    if(selection->getSwitch("ecalIso") && 
       !selection->passCut("ecalIso", m_ecalIso)) return false; 
    m_electronCounter.IncrVar("ecalIso");
    
    if(selection->getSwitch("trackerIso") && 
       !selection->passCut("trackerIso", m_trackerIso)) return false; 
    m_electronCounter.IncrVar("trackerIso");
    
    if(selection->getSwitch("hcalIso") && 
       !selection->passCut("hcalIso", m_hcalIso)) return false; 
    m_electronCounter.IncrVar("hcalIso");
  }    

  if(selection->getSwitch("combIso")) {
    if(!selection->passCut("combIso", m_combIso)) return false; 
    m_electronCounter.IncrVar("combIso");
  }

  if(selection->getSwitch("combPFIso")) {
    if(!selection->passCut("combPFIso", m_combPFIso)) return false; 
    m_electronCounter.IncrVar("combPFIso");
  }

  m_electronCounter.IncrVar("finalCustomEleIDOnlyIso");

  return true;
}

bool CutBasedEleIDSelector::outputConv()
{
  if ( m_fiducialflag == -1 ) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CutBasedEleIDSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return false;
  }

  if ( !m_classDep ) {
    cout << "CutBasedEleIDSelector::output() method is intended for class based ele ID, while you configured for class independent ID. Use outputNoClass() instead." << endl;
    return false;
  }

  int GsfClass = m_electronClass;
  int offset=0;
  Utils anaUtils;
  bool eeElectron = anaUtils.fiducialFlagECAL(m_fiducialflag, isEE);
  if(eeElectron) {
    offset=4;
  }

  Selection *selection;
  if(GsfClass==GOLDEN) 
    selection=m_EgammaCutIDSelection[offset];
  else if(GsfClass==BIGBREM)
    selection=m_EgammaCutIDSelection[offset+1];
  else if(GsfClass==NARROW)
    selection=m_EgammaCutIDSelection[offset+2];
  else if(GsfClass==SHOWERING) {
    selection=m_EgammaCutIDSelection[offset+3];
  }

  m_electronCounter.IncrVar("electronsOnlyConv");

  if(selection->getSwitch("missHits") && 
     !selection->passCut("missHits", m_missingHits)) return false; 
  m_electronCounter.IncrVar("missHits");

  if(selection->getSwitch("distConv") && selection->getSwitch("dcotConv") &&
     (!selection->passCut("distConv", m_distConv) && !selection->passCut("dcotConv", m_dcotConv)) ) return false; 
  m_electronCounter.IncrVar("vetoConv");

  if(selection->getSwitch("matchedConv") && 
     m_matchedConv) return false; 
  m_electronCounter.IncrVar("matchedConv");

  m_electronCounter.IncrVar("finalCustomEleIDOnlyConv");

  return true;
}

bool CutBasedEleIDSelector::outputNoClass() 
{
  if ( m_fiducialflag == -1 ) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CutBasedEleIDSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return false;
  }

  if ( m_classDep ) {
    cout << "CutBasedEleIDSelector::output() method is intended for class IN-dependent ele ID, while you configured for class based ID. Use output() instead." << endl;
    return false;
  }

  m_electronCounter.IncrVar("electrons");

  Utils anaUtils;
  bool isPFlowOnly = !anaUtils.electronRecoType(m_recoflag, isEcalDriven);
  if ( isPFlowOnly && !m_applyIDOnPFlow ) return true;

  if (!outputNoClassEleId())
    return false;
  if (!outputNoClassIso())
    return false;
  if (!outputNoClassConv())
    return false;
  m_electronCounter.IncrVar("finalCustomEleID");

  return true;  

}

bool CutBasedEleIDSelector::outputNoClassEleId() {

  if ( m_fiducialflag == -1 ) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CutBasedEleIDSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return false;
  }

  if ( m_classDep ) {
    cout << "CutBasedEleIDSelector::output() method is intended for class IN-dependent ele ID, while you configured for class based ID. Use output() instead." << endl;
    return false;
  }

  m_electronCounter.IncrVar("electronsOnlyID");

  Utils anaUtils;
  bool isPFlowOnly = !anaUtils.electronRecoType(m_recoflag, isEcalDriven);
  if ( isPFlowOnly && !m_applyIDOnPFlow ) return true;

  bool eleInGap = anaUtils.fiducialFlagECAL(m_fiducialflag, isEBEEGap);

  //  if ( eleInGap ) return false;

  int offset=0;
  bool eeElectron = anaUtils.fiducialFlagECAL(m_fiducialflag, isEE);
  if(eeElectron) offset=1;

  Selection *selection = m_EgammaCutIDSelection[offset];

  if(selection->getSwitch("egammaCutBased") && 
     m_egammaCutBased ) m_electronCounter.IncrVar("egammaCutBased");

  if(selection->getSwitch("likelihood0") && 
     ( m_NBrem == 0 && !selection->passCut("likelihood0", m_Likelihood ) ) 
     ) return false;
  if(selection->getSwitch("likelihood1") && 
     ( m_NBrem > 0 && !selection->passCut("likelihood1", m_Likelihood ) ) 
     ) return false;
  m_electronCounter.IncrVar("likelihood");
  
  if(selection->getSwitch("hOverE") && 
     !selection->passCut("hOverE", fabs(m_HOverE) )) return false;  
  m_electronCounter.IncrVar("hOverE");
  
  if(selection->getSwitch("s9s25") && 
     !selection->passCut("s9s25", m_S9S25)) return false; 
  m_electronCounter.IncrVar("s9s25");
  
  if(selection->getSwitch("deta") && 
     !selection->passCut("deta", fabs(m_DEta) )) return false; 
  m_electronCounter.IncrVar("deta");

  if(selection->getSwitch("dphiIn") && 
     !selection->passCut("dphiIn", fabs(m_DPhiIn) )) return false; 
  m_electronCounter.IncrVar("dphiIn");

  if(selection->getSwitch("dphiOut") && 
     !selection->passCut("dphiOut", fabs(m_DPhiOut) )) return false; 
  m_electronCounter.IncrVar("dphiOut");
    
  if(selection->getSwitch("invEminusInvP") && 
     !selection->passCut("invEminusInvP", m_InvEminusInvP)) return false; 
  m_electronCounter.IncrVar("invEminusInvP");
    
  if(selection->getSwitch("bremFraction") && 
     !selection->passCut("bremFraction", m_BremFraction)) return false; 
  m_electronCounter.IncrVar("bremFraction");

  if(selection->getSwitch("sigmaEtaEta") && 
     !selection->passCut("sigmaEtaEta", m_SigmaEtaEta)) return false; 
  m_electronCounter.IncrVar("sigmaEtaEta");

  if(selection->getSwitch("sigmaPhiPhi") && 
     !selection->passCut("sigmaPhiPhi", m_SigmaPhiPhi)) return false; 
  m_electronCounter.IncrVar("sigmaPhiPhi");

  if(selection->getSwitch("eOverPout") && 
     !selection->passCut("eOverPout", m_EOverPout)) return false; 
  m_electronCounter.IncrVar("eOverPout");

  if(selection->getSwitch("eOverPin") && 
     !selection->passCut("eOverPin", m_EOverPin)) return false; 
  m_electronCounter.IncrVar("eOverPin");

  if(m_doEcalCleaning && !m_cleaner->output()) return false;
  m_electronCounter.IncrVar("ecalCleaning");

  m_electronCounter.IncrVar("finalCustomEleIDOnlyID");

  return true;

}

bool CutBasedEleIDSelector::outputNoClassIso() 
{
  if ( m_fiducialflag == -1 ) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CutBasedEleIDSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return false;
  }

  if ( m_classDep ) {
    cout << "CutBasedEleIDSelector::output() method is intended for class IN-dependent ele ID, while you configured for class based ID. Use output() instead." << endl;
    return false;
  }

  int offset=0;
  Utils anaUtils;
  bool eeElectron = anaUtils.fiducialFlagECAL(m_fiducialflag, isEE);
  if(eeElectron) offset=1;

  Selection *selection = m_EgammaCutIDSelection[offset];

  m_electronCounter.IncrVar("electronsOnlyIso");

  if(!selection->getSwitch("combIso") && !selection->getSwitch("combPFIso")) {
    if(selection->getSwitch("ecalIso") && 
       !selection->passCut("ecalIso", m_ecalIso)) return false; 
    m_electronCounter.IncrVar("ecalIso");
    
    if(selection->getSwitch("trackerIso") && 
       !selection->passCut("trackerIso", m_trackerIso)) return false; 
    m_electronCounter.IncrVar("trackerIso");
    
    if(selection->getSwitch("hcalIso") && 
       !selection->passCut("hcalIso", m_hcalIso)) return false; 
    m_electronCounter.IncrVar("hcalIso");
  }

  if(selection->getSwitch("combIso")) {
    if( !selection->passCut("combIso", m_combIso)) return false; 
    m_electronCounter.IncrVar("combIso");
  }

  if(selection->getSwitch("combPFIso")) {
    if( !selection->passCut("combPFIso", m_combPFIso)) return false; 
    m_electronCounter.IncrVar("combPFIso");
  }

  m_electronCounter.IncrVar("finalCustomEleIDOnlyIso");

  return true;
}

bool CutBasedEleIDSelector::outputNoClassConv()
{
  if ( m_fiducialflag == -1 ) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CutBasedEleIDSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return false;
  }

  if ( m_classDep ) {
    cout << "CutBasedEleIDSelector::output() method is intended for class IN-dependent ele ID, while you configured for class based ID. Use output() instead." << endl;
    return false;
  }

  int offset=0;
  Utils anaUtils;
  bool eeElectron = anaUtils.fiducialFlagECAL(m_fiducialflag, isEE);
  if(eeElectron) offset=1;

  Selection *selection = m_EgammaCutIDSelection[offset];

  m_electronCounter.IncrVar("electronsOnlyConv");

  if(selection->getSwitch("missHits") && 
     !selection->passCut("missHits", m_missingHits)) return false; 
  m_electronCounter.IncrVar("missHits");

  if(selection->getSwitch("distConv") && selection->getSwitch("dcotConv") &&
     (!selection->passCut("distConv", m_distConv) && !selection->passCut("dcotConv", m_dcotConv)) ) return false; 
  m_electronCounter.IncrVar("vetoConv");

  if(selection->getSwitch("matchedConv") && 
     m_matchedConv) return false; 
  m_electronCounter.IncrVar("matchedConv");
  
  m_electronCounter.IncrVar("finalCustomEleIDOnlyConv");

  return true;
}

void CutBasedEleIDSelector::displayEfficiencies() {

  m_electronCounter.Draw();

  std::cout << "++++++ ELEID PART +++++++" << std::endl;
  m_electronCounter.Draw("finalCustomEleIDOnlyID","electronsOnlyID");
  std::cout << "====== CUTS EFFICIENCY ======" << std::endl;
  m_electronCounter.Draw("hOverE","electronsOnlyID");
  m_electronCounter.Draw("s9s25","hOverE");
  m_electronCounter.Draw("deta","s9s25");
  m_electronCounter.Draw("dphiIn","deta");
  m_electronCounter.Draw("dphiOut","dphiIn");
  m_electronCounter.Draw("invEminusInvP", "dphiOut");
  m_electronCounter.Draw("bremFraction", "invEminusInvP");
  m_electronCounter.Draw("sigmaEtaEta","bremFraction");
  m_electronCounter.Draw("sigmaPhiPhi","sigmaEtaEta");
  m_electronCounter.Draw("eOverPout", "sigmaPhiPhi");
  m_electronCounter.Draw("eOverPin", "eOverPout");
  m_electronCounter.Draw("ecalCleaning", "eOverPin");

  std::cout << "++++++ ISO PART +++++++" << std::endl;
  m_electronCounter.Draw("finalCustomEleIDOnlyIso","electronsOnlyIso");
  std::cout << "====== CUTS EFFICIENCY ======" << std::endl;
  m_electronCounter.Draw("ecalIso", "electronsOnlyIso");
  m_electronCounter.Draw("trackerIso", "ecalIso");
  m_electronCounter.Draw("hcalIso", "trackerIso");
  m_electronCounter.Draw("combIso", "electronsOnlyIso");
  m_electronCounter.Draw("combPFIso", "combIso");

  std::cout << "++++++ CONV REJ PART +++++++" << std::endl;
  m_electronCounter.Draw("finalCustomEleIDOnlyConv","electronsOnlyConv");
  std::cout << "====== CUTS EFFICIENCY ======" << std::endl;
  m_electronCounter.Draw("missHits", "electronsOnlyConv");
  m_electronCounter.Draw("vetoConv", "missHits");
  m_electronCounter.Draw("matchedConv", "vetoConv");

  std::cout << "++++++ FINAL EFFICIENCY +++++++" << std::endl;
  m_electronCounter.Draw("finalCustomEleID","electrons");
  m_electronCounter.Draw("likelihood","electrons");
  m_electronCounter.Draw("egammaCutBased", "electrons");

  if(m_doEcalCleaning) m_cleaner->displayEfficiencies();
}
