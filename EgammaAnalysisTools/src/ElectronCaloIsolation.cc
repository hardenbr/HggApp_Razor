#include <iostream>
#include "EgammaAnalysisTools/include/ElectronCaloIsolation.hh"

ElectronCaloIsolation::ElectronCaloIsolation (){}

ElectronCaloIsolation::ElectronCaloIsolation (TLorentzVector electronSuperCluster) {
  m_electronSuperCluster = electronSuperCluster;

  //! default configuration
  m_extRadius = 0.40; // deltaR
  m_intRadius = 0.00; // deltaR
  m_relative = true;
  m_electronAtVertex.SetPxPyPzE(0,0,0,0);
  
  m_EtHcal = 0.0;
  m_EtEcal = 0.0;

}

void ElectronCaloIsolation::setCaloTowers(std::vector<TVector3> position, 
					  std::vector<float> emEnergy, 
					  std::vector<float> hadEnergy) { 

  m_towersPosition = position;
  m_emEnergy = emEnergy;
  m_hadEnergy = hadEnergy;

}

void ElectronCaloIsolation::getEtTowers() {

  float tmpSumHadEt = 0.0;
  float tmpSumEmEt = 0.0;

  std::vector<TVector3>::const_iterator towerPosition;
  std::vector<float>::const_iterator theEmEnergy = m_emEnergy.begin();
  std::vector<float>::const_iterator theHadEnergy = m_hadEnergy.begin();
  for(towerPosition=m_towersPosition.begin(); towerPosition!=m_towersPosition.end(); 
      towerPosition++) { 
    
    float hadEt = *theHadEnergy;
    float emEt = *theEmEnergy;
    
    TVector3 scPosition = m_electronSuperCluster.Vect();
    double dr = scPosition.DeltaR(*towerPosition);

    if ( fabs(dr) < m_extRadius && 
	 fabs(dr) > m_intRadius ){ 
      
      tmpSumHadEt += hadEt;
      tmpSumEmEt += emEt;

    }
    
    theHadEnergy++;
    theEmEnergy++;

  }

  m_EtHcal = tmpSumHadEt;
  m_EtEcal = tmpSumEmEt;

}

float ElectronCaloIsolation::getSumEtHcal() {

  getEtTowers();

  float xval=0;
  float ele_ptAtVtx = m_electronAtVertex.Pt();

  if( m_relative && ele_ptAtVtx == 0 ) {
    std::cout << "relative calo isolation required. But electron momentum at vertex not set."
	      << "\nPlease use setElectronMomentumAtVtx(TLorentzVector electronAtVertex)!"
	      << "\nreturning sumEt=10000" << std::endl;
    return 10000.;
  }
  else if( m_relative && ele_ptAtVtx!=0 ) {
    double elePt = m_electronAtVertex.Pt();
    xval = m_EtHcal / elePt;
  }
  else xval = m_EtHcal;

  return xval;

}

float ElectronCaloIsolation::getSumEtEcal() {

  getEtTowers();

  float xval=0;
  float ele_ptAtVtx = m_electronAtVertex.Pt();

  if( m_relative && ele_ptAtVtx == 0 ) {
    std::cout << "relative calo isolation required. But electron momentum at vertex not set."
	      << "\nPlease use setElectronMomentumAtVtx(TLorentzVector electronAtVertex)!"
	      << "\nreturning sumEt=10000" << std::endl;
    return 10000.;
  }
  else if( m_relative && ele_ptAtVtx!=0 ) {
    double elePt = m_electronAtVertex.Pt();
    xval = m_EtEcal / elePt;
  }
  else xval = m_EtEcal;

  return xval;

}

