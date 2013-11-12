#include "EgammaAnalysisTools/include/EcalCleaner.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/EcalRecHitBits.h"
#include "CommonTools/include/Utils.hh"
#include <iostream>
#include <math.h>

using namespace bits;

EcalCleaner::EcalCleaner() {}

EcalCleaner::~EcalCleaner() {}

void EcalCleaner::Configure(const char *configDir) 
{

  std::string fileSwitches = std::string(configDir) + "/ecalCleaningSwitches.txt";
  std::string cuts = std::string(configDir) + "/ecalCleaningCuts.txt";
  
  m_selection = new Selection(cuts, fileSwitches);

  m_selection->addCut("spikeFraction");
  m_selection->addCut("seedTime");
  m_selection->addSwitch("removeWeirdSeed");
  m_selection->summary();

  m_electronCounter.SetTitle("SINGLE CLUSTER COUNTER");
  m_electronCounter.AddVar("barrelElectrons");
  m_electronCounter.AddVar("spikeFraction");
  m_electronCounter.AddVar("seedTime");
  m_electronCounter.AddVar("removeWeirdSeed");
  m_electronCounter.AddVar("finalEcalCleaning");

}

bool EcalCleaner::output() {

  m_electronCounter.IncrVar("electronsOnlyID");

  Utils anaUtils;

  bool eleInBarrel = anaUtils.fiducialFlagECAL(m_fiducialFlag, isEB);

  if ( !eleInBarrel ) return true;
  m_electronCounter.IncrVar("barrelElectrons");

  if(m_selection->getSwitch("spikeFraction") && 
     !m_selection->passCut("spikeFraction", 1.0-m_e4SwissCross/m_e1) ) return false;  
  m_electronCounter.IncrVar("spikeFraction");
  
  if(m_selection->getSwitch("seedTime") && 
     !m_selection->passCut("seedTime", m_seedTime) ) return false;  
  m_electronCounter.IncrVar("seedTime");

  if(m_selection->getSwitch("removeWeirdSeed") && 
     m_seedFlag != EcalRecHitBits::kGood ) return false;  
  m_electronCounter.IncrVar("removeWeirdSeed");

  m_electronCounter.IncrVar("finalEcalCleaning");

  return true;

}

void EcalCleaner::displayEfficiencies() {

  std::cout << "====== ECAL CLEANING CUTS EFFICIENCY ======" << std::endl;
  m_electronCounter.Draw();
  m_electronCounter.Draw("spikeFraction","barrelElectrons");
  m_electronCounter.Draw("seedTime","spikeFraction");
  m_electronCounter.Draw("removeWeirdSeed","seedTime");
  m_electronCounter.Draw("finalEcalCleaning","barrelElectrons");

}
