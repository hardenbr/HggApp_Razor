#include <algorithm>
#include "EgammaAnalysisTools/include/ElectronBestCandidateSelector.hh" 

ElectronBestCandidateSelector::ElectronBestCandidateSelector( std::vector<ElectronQualityData> electrons ) {

  m_electrons = electrons;

}

std::pair<int,int> ElectronBestCandidateSelector::bestByPt() {

  if(m_electrons.size()==0) return std::make_pair(-1,-1);
  
  if(m_electrons.size()==1) return std::make_pair( m_electrons.begin()->index, -1);

  sort( m_electrons.begin(), m_electrons.end(), PtSorting() );
  int first = m_electrons.begin()->index;
  int second = (m_electrons.begin()+1)->index;
  std::pair<int,int> bestPair = std::make_pair(first,second);
  return bestPair;

}


std::pair<int,int> ElectronBestCandidateSelector::bestBySCenergy() {

  if(m_electrons.size()==0) return std::make_pair(-1,-1);
  
  if(m_electrons.size()==1) return std::make_pair( m_electrons.begin()->index, -1);

  sort( m_electrons.begin(), m_electrons.end(), SCenergySorting() );
  int first = m_electrons.begin()->index;
  int second = (m_electrons.begin()+1)->index;
  std::pair<int,int> bestPair = std::make_pair(first,second);
  return bestPair;

}


std::pair<int,int> ElectronBestCandidateSelector::bestByTrackerIsolation() {

  if(m_electrons.size()==0) return std::make_pair(-1,-1);
  
  if(m_electrons.size()==1) return std::make_pair( m_electrons.begin()->index, -1);

  sort( m_electrons.begin(), m_electrons.end(), TkIsolationSorting() );
  int first = m_electrons.begin()->index;
  int second = (m_electrons.begin()+1)->index;
  std::pair<int,int> bestPair = std::make_pair(first,second);
  return bestPair;

}


std::pair<int,int> ElectronBestCandidateSelector::bestByEcalIsolation() {

  if(m_electrons.size()==0) return std::make_pair(-1,-1);
  
  if(m_electrons.size()==1) return std::make_pair( m_electrons.begin()->index, -1);

  sort( m_electrons.begin(), m_electrons.end(), EcalIsolationSorting() );
  int first = m_electrons.begin()->index;
  int second = (m_electrons.begin()+1)->index;
  std::pair<int,int> bestPair = std::make_pair(first,second);
  return bestPair;

}


std::pair<int,int> ElectronBestCandidateSelector::bestByHcalIsolation() {

  if(m_electrons.size()==0) return std::make_pair(-1,-1);
  
  if(m_electrons.size()==1) return std::make_pair( m_electrons.begin()->index, -1);

  sort( m_electrons.begin(), m_electrons.end(), HcalIsolationSorting() );
  int first = m_electrons.begin()->index;
  int second = (m_electrons.begin()+1)->index;
  std::pair<int,int> bestPair = std::make_pair(first,second);
  return bestPair;

}


std::pair<int,int> ElectronBestCandidateSelector::bestByElectronIdLH() {

  if(m_electrons.size()==0) return std::make_pair(-1,-1);
  
  if(m_electrons.size()==1) return std::make_pair( m_electrons.begin()->index, -1);

  sort( m_electrons.begin(), m_electrons.end(), LHSorting() );
  int first = m_electrons.begin()->index;
  int second = (m_electrons.begin()+1)->index;
  std::pair<int,int> bestPair = std::make_pair(first,second);
  return bestPair;

}

