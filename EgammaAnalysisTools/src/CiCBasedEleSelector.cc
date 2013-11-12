#include "EgammaAnalysisTools/include/CiCBasedEleSelector.hh"
#include "EgammaAnalysisTools/include/eIDCiCCuts.h"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/Utils.hh"
#include <iostream>
#include <math.h>
#include <string.h>

using namespace bits;

CiCBasedEleSelector::CiCBasedEleSelector() {

  // if not set, assign the flag to ECAL driven
  m_recoflag = 2;

  m_doEcalCleaning = false;

  // fiducial flag HAS to be set, otherwise the selector doesn't know which slection apply (EB or EE)
  m_fiducialflag = -1;

  m_useClass = false;
}

CiCBasedEleSelector::~CiCBasedEleSelector() {}

void CiCBasedEleSelector::Configure(std::string type, bool useEtBins, bool specialCategories, int version)
{

  if (type == "CiCVeryLoose")
    m_eIDLevel = 0;
  else if (type == "CiCLoose")
    m_eIDLevel = 1;
  else if (type == "CiCMedium")
    m_eIDLevel = 2;
  else if (type == "CiCTight")
    m_eIDLevel = 3;
  else if (type == "CiCSuperTight")
    m_eIDLevel = 4;
  else if (type == "CiCHyperTight")
    m_eIDLevel = 5;
  else if (type == "CiCHyperTight2")
    m_eIDLevel = 6;
  else if (type == "CiCHyperTight3")
    m_eIDLevel = 7;
  else if (type == "CiCHyperTight4")
    m_eIDLevel = 8;

  m_version = version;
  m_useEtBins = useEtBins;
  m_specialCategories = specialCategories;

  m_electronCounter.SetTitle("SINGLE ELECTRON COUNTER");
  m_electronCounter.AddVar("electrons");
  m_electronCounter.AddVar("electronsOnlyID");
  m_electronCounter.AddVar("electronsOnlyIso");
  m_electronCounter.AddVar("electronsOnlyConv");
  m_electronCounter.AddVar("finalCustomEleID");
  m_electronCounter.AddVar("finalCustomEleIDOnlyID");
  m_electronCounter.AddVar("finalCustomEleIDOnlyIso");
  m_electronCounter.AddVar("finalCustomEleIDOnlyConv");
}

void CiCBasedEleSelector::ConfigureEcalCleaner(const char *configDir) {

  m_doEcalCleaning = true;

  m_cleaner = new EcalCleaner();
  m_cleaner->Configure(configDir);

}


bool CiCBasedEleSelector::output() {
  ElectronClassification();

  m_electronCounter.IncrVar("electrons");

  if (!outputEleId())
    return false;
  if (!outputIso())
    return false;
  if (!outputConv())
    return false;
  m_electronCounter.IncrVar("finalCustomEleID");

  return true;  

}

void CiCBasedEleSelector::ElectronClassification() 
{
  if ( m_fiducialflag == -1) {
    cout << "CONFIGURATION ERROR! Fiducial flag not set. Use the method CiCBasedEleSelector::SetEcalFiducialRegion(int word) to do it. Returning always false eleID." << endl;
    return;
  }

  Utils anaUtils;
  bool isPFlowOnly = !anaUtils.electronRecoType(m_recoflag, isEcalDriven);
  
  bool ebElectron = anaUtils.fiducialFlagECAL(m_fiducialflag, isEB);
  bool eeElectron = anaUtils.fiducialFlagECAL(m_fiducialflag, isEE);
  
  float fbrem = m_BremFraction;
  float eopin = m_EOverPin;

  if (ebElectron) {
    if (m_BremFraction>=0.12 && m_EOverPin >0.9 && m_EOverPin < 1.2)    // bremming barrel
      m_cat = 0;
    else if (m_specialCategories && ((fabs(m_SCEta) >  .445   && fabs(m_SCEta) <  .45  ) ||      // barrel crack electrons
                                   (fabs(m_SCEta) >  .79    && fabs(m_SCEta) <  .81  ) ||
                                   (fabs(m_SCEta) > 1.137   && fabs(m_SCEta) < 1.157 ) ||
                                   (fabs(m_SCEta) > 1.47285 && fabs(m_SCEta) < 1.4744)))
      m_cat = 6;                                       
    else if (isPFlowOnly)              // tracker driven electrons
      m_cat = 8;
    else if (m_BremFraction < 0.12)                           // low brem barrel 
      m_cat = 1;
    else                                             // bad track in barrel
      m_cat = 2;
  } else {
    if (m_BremFraction>=0.2 && m_EOverPin >0.82 && m_EOverPin < 1.22)   // bremming endcap
      m_cat = 3;
    else if (m_specialCategories && (fabs(m_SCEta) > 1.5 && fabs(m_SCEta) < 1.58))                // endcap crack electrons
      m_cat = 7;                                       
    else if (isPFlowOnly)              // tracker driven electrons
      m_cat = 8;
    else if (m_BremFraction < 0.2)                            // low brem endcap 
      m_cat = 4;
    else                                             // bad track in barrel
      m_cat = 5;
  }

  if (m_useEtBins)
    {
      if ( m_SCEt < 20.) 
	m_ptBin = 2;
      else if ( m_SCEt > 30.)
	m_ptBin = 0;
      else
	m_ptBin = 1;
    }  
  
  m_useClass=true;

}

bool CiCBasedEleSelector::reset()
{
  m_useClass=false;
  m_ptBin=-1;
  m_cat=-1;
  m_SCEt=-1;
  m_SCEta=-1;
  m_ESeedOverPin=-1;
  m_HOverE=-1;
  m_DEta=-1;
  m_DPhiIn=-1;
  m_BremFraction=-1;
  m_SigmaEtaEta=-1;
  m_EOverPin=-1;
  m_ecalIso=-1;
  m_trackerIso=-1;
  m_hcalIso=-1;
  m_distConv=-1;
  m_dcotConv=-1;
  m_missingHits=-1;
  m_fiducialflag=-1;
  m_recoflag=-1;
}

bool CiCBasedEleSelector::outputEleId() {

  ElectronClassification();


  m_electronCounter.IncrVar("electronsOnlyID");

  float eseedopincor = m_ESeedOverPin + m_BremFraction;
  if(m_BremFraction<0) 
    eseedopincor = m_ESeedOverPin; 

  if (m_BremFraction < -2) {
    return false;
  }
  

   if (
       m_version == 1 && 
       (m_HOverE < cuthoe[m_ptBin][m_eIDLevel][m_cat]) &&
       (m_SigmaEtaEta < cutsee[m_ptBin][m_eIDLevel][m_cat]) &&
       (fabs(m_DPhiIn) < cutdphiin[m_ptBin][m_eIDLevel][m_cat])  &&
       (fabs(m_DEta) < cutdetain[m_ptBin][m_eIDLevel][m_cat]) &&
       (eseedopincor > cuteseedopcor[m_ptBin][m_eIDLevel][m_cat]) &&
       (m_SCEt > cutet[m_ptBin][m_eIDLevel][m_cat])
       ) 
     {
      m_electronCounter.IncrVar("finalCustomEleIDOnlyID");  
      return true;
     } 
   else
     if (
	 m_version == 2 && 
	 compute_eid_cut(m_HOverE, m_SCEt, cuthoelv03_old[0][m_eIDLevel][m_cat], cuthoev03_old[0][m_eIDLevel][m_cat], false) &&
	 compute_eid_cut(m_SigmaEtaEta, m_SCEt, cutseelv03_old[0][m_eIDLevel][m_cat], cutseev03_old[0][m_eIDLevel][m_cat], false) &&
	 compute_eid_cut(fabs(m_DPhiIn), m_SCEt, cutdphiinlv03_old[0][m_eIDLevel][m_cat], cutdphiinv03_old[0][m_eIDLevel][m_cat], false) &&
	 compute_eid_cut(fabs(m_DEta), m_SCEt, cutdetainlv03_old[0][m_eIDLevel][m_cat], cutdetainv03_old[0][m_eIDLevel][m_cat], false) &&
	 eseedopincor > cuteseedopcorv03_old[0][m_eIDLevel][m_cat] 
	 )
       {
	 m_electronCounter.IncrVar("finalCustomEleIDOnlyID");  
	 return true;
       } 
   else
     if (
	 //Fixed isolation version
	 m_version == 3 && 
	 compute_eid_cut(m_HOverE, m_SCEt, cuthoelv04[0][m_eIDLevel][m_cat], cuthoev04[0][m_eIDLevel][m_cat], false) &&
	 compute_eid_cut(m_SigmaEtaEta, m_SCEt, cutseelv04[0][m_eIDLevel][m_cat], cutseev04[0][m_eIDLevel][m_cat], false) &&
	 compute_eid_cut(fabs(m_DPhiIn), m_SCEt, cutdphiinlv04[0][m_eIDLevel][m_cat], cutdphiinv04[0][m_eIDLevel][m_cat], false) &&
	 compute_eid_cut(fabs(m_DEta), m_SCEt, cutdetainlv04[0][m_eIDLevel][m_cat], cutdetainv04[0][m_eIDLevel][m_cat], false) &&
	 eseedopincor > cuteseedopcorv04[0][m_eIDLevel][m_cat] 
	 )
       {
	 m_electronCounter.IncrVar("finalCustomEleIDOnlyID");  
	 return true;
       } 


  
  return false;
}


bool CiCBasedEleSelector::outputIso() 
{

  ElectronClassification();

  m_electronCounter.IncrVar("electronsOnlyIso");
  
  float iso_sum = m_trackerIso + m_ecalIso + m_hcalIso;   

  if(fabs(m_SCEta)>1.5) 
    iso_sum += (fabs(m_SCEta)-1.5)*1.09;

  float iso_sumoet = iso_sum*(40./m_SCEt);
  float iso_sum_scaled = iso_sum*pow(40./m_SCEt, 2); 

  if ( m_version == 1 &&
       (iso_sum < cutiso_sum[m_ptBin][m_eIDLevel][m_cat]) &&
       (iso_sum_scaled < cutiso_sumoet[m_ptBin][m_eIDLevel][m_cat])
       )
    {
      m_electronCounter.IncrVar("finalCustomEleIDOnlyIso");  
      return true;
    }
  else if ( m_version == 2 && 
	    compute_eid_cut(iso_sumoet, m_SCEt, cutiso_sumoetlv03_old[0][m_eIDLevel][m_cat], cutiso_sumoetv03_old[0][m_eIDLevel][m_cat], false) &&
	    (iso_sum < cutiso_sumv03_old[0][m_eIDLevel][m_cat])
	    )
    {
      m_electronCounter.IncrVar("finalCustomEleIDOnlyIso");  
      return true;
    }
  else if ( m_version == 3 && 
	    compute_eid_cut(iso_sumoet, m_SCEt, cutiso_sumoetlv04[0][m_eIDLevel][m_cat], cutiso_sumoetv04[0][m_eIDLevel][m_cat], false) &&
	    (iso_sum < cutiso_sumv04[0][m_eIDLevel][m_cat])
	    )
    {
      m_electronCounter.IncrVar("finalCustomEleIDOnlyIso");  
      return true;
    }
  
  return false;
}

bool CiCBasedEleSelector::outputConv()
{
  ElectronClassification();
  
  m_electronCounter.IncrVar("electronsOnlyConv");
  
  float dcotdist = ((0.04 - std::max(fabs(m_distConv), fabs(m_dcotConv))) > 0?(0.04 - std::max(fabs(m_distConv), fabs(m_dcotConv))):0);
  
  if ( m_version == 1 &&
       (float(m_missingHits) < cutfmishits[m_ptBin][m_eIDLevel][m_cat]) && (dcotdist < cutdcotdist[m_ptBin][m_eIDLevel][m_cat])
       )
    {
      m_electronCounter.IncrVar("finalCustomEleIDOnlyConv");  
      return true;
    }
  else
    if ( m_version == 2 &&
	 float(m_missingHits) < cutfmishitsv03_old[0][m_eIDLevel][m_cat]
	 )
      {
	m_electronCounter.IncrVar("finalCustomEleIDOnlyConv");  
	return true;
      } 
  else
    if ( m_version == 3 &&
	 float(m_missingHits) < cutfmishitsv04[0][m_eIDLevel][m_cat]
	 )
      {
	m_electronCounter.IncrVar("finalCustomEleIDOnlyConv");  
	return true;
      } 
  
  return false;
}

void CiCBasedEleSelector::displayEfficiencies() {
  m_electronCounter.Draw();
  std::cout << "++++++ ELEID PART +++++++" << std::endl;
  m_electronCounter.Draw("finalCustomEleIDOnlyID","electronsOnlyID");
  std::cout << "++++++ ISO PART +++++++" << std::endl;
  m_electronCounter.Draw("finalCustomEleIDOnlyIso","electronsOnlyIso");
  std::cout << "++++++ CONV REJ PART +++++++" << std::endl;
  m_electronCounter.Draw("finalCustomEleIDOnlyConv","electronsOnlyConv");
  std::cout << "++++++ FINAL EFFICIENCY +++++++" << std::endl;
  m_electronCounter.Draw("finalCustomEleID","electrons");
  if(m_doEcalCleaning) m_cleaner->displayEfficiencies();
}

// Linearization of the cut
bool CiCBasedEleSelector::compute_eid_cut(float x, float et, float cut_min, float cut_max, bool gtn) {

  float et_min = 10;
  float et_max = 40;

  bool accept = false;
  float cut = cut_max;    //  the cut at et=40 GeV

  if(et < et_max) {
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  } 

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  return accept;
}
