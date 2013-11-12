#include "EgammaAnalysisTools/include/ElectronTrackerIsolation.hh"

ElectronTrackerIsolation::ElectronTrackerIsolation (){}

ElectronTrackerIsolation::ElectronTrackerIsolation (TLorentzVector electron, float lipEle) {
  m_electron = electron;
  m_lipEle = lipEle;

  //! default configuration
  m_extRadius = 0.40; // deltaR
  m_intRadius = 0.02; // deltaR
  m_maxDeltaLip = 0.2; // cm
  m_minPtTrack = 1.5; // GeV
  m_relative = true;
}

ElectronTrackerIsolation::~ElectronTrackerIsolation (){}

void ElectronTrackerIsolation::setTracks( std::vector<TLorentzVector> tracks, std::vector<float> lips) {
  m_tracks = tracks;
  m_lips = lips;
}

float ElectronTrackerIsolation::getSumPtTracks() {

  float dummyPt = 0 ;

  float ele_pt  = m_electron.Pt();

  std::vector<float>::const_iterator trackLip = m_lips.begin();
  std::vector<TLorentzVector>::const_iterator track;
  for(track=m_tracks.begin(); track!=m_tracks.end(); track++) { 
    
    float track_pt = track->Pt();

    // only tracks from the same vertex as the electron
    float track_lip = *trackLip;

    if ( fabs(track_lip - m_lipEle) > m_maxDeltaLip ) continue;

    // usually low pt tracks are fakes
    if ( track_pt < m_minPtTrack ) continue;

    double dr = m_electron.DeltaR(*track);
    if ( fabs(dr) < m_extRadius && fabs(dr) > m_intRadius ){ 
      if ( m_relative ) {
	dummyPt += track_pt/ele_pt;
      }
      else {
	dummyPt += track_pt;
      }
    } 
    
    trackLip++;
  } //end loop over tracks		       
  
  return dummyPt;
}


float ElectronTrackerIsolation::getSumPtSquaredTracks() {

  float dummyPt = 0 ;

  float ele_pt  = m_electron.Pt();

  std::vector<float>::const_iterator trackLip = m_lips.begin();
  std::vector<TLorentzVector>::const_iterator track;
  for(track=m_tracks.begin(); track!=m_tracks.end(); track++) { 
    
    float track_pt = track->Pt();

    // only tracks from the same vertex as the electron
    float track_lip = *trackLip;

    if ( fabs(track_lip - m_lipEle) > m_maxDeltaLip ) continue;

    // usually low pt tracks are fakes
    if ( track_pt < m_minPtTrack ) continue;

    double dr = m_electron.DeltaR(*track);
    if ( fabs(dr) < m_extRadius && fabs(dr) > m_intRadius ){ 
      if ( m_relative ) {
	dummyPt += pow(track_pt/ele_pt,2);
      }
      else {
	dummyPt += pow(track_pt,2);
      }
    } 
    
    trackLip++;
  } //end loop over tracks		       
  
  return dummyPt;
}

