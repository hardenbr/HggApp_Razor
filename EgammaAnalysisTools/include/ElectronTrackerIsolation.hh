#ifndef ELECTRONTRACKERISOLATION_H
#define ELECTRONTRACKERISOLATION_H

#include <vector>
#include "TLorentzVector.h"

class ElectronTrackerIsolation{
public:
  
  //! constructors
  ElectronTrackerIsolation();
  ElectronTrackerIsolation( TLorentzVector electron, float lipEle);

  //! destructor 
  ~ElectronTrackerIsolation();

  //! tracks of the events + lips
  void setTracks( std::vector<TLorentzVector> tracks, std::vector<float> lips);

  //! internal and external cones
  void setExtRadius (float extRadius) { m_extRadius = extRadius; } 
  void setIntRadius (float intRadius) { m_intRadius = intRadius; }
  //! set deltaz consistency of the vertex
  void setMaxLip(float maxDeltaLip) { m_maxDeltaLip = maxDeltaLip; }
  //! set the min pt of the track to be considered
  void setMinPtTrack(float minPtTrack) { m_minPtTrack = minPtTrack; }
  //! if true, returns sumPt / ptEle
  void setRelative(bool relative) { m_relative = relative; }

  //! returns sumPtTracks or sumPtTracks/ptEle according to m_relative
  float getSumPtTracks();
  //! returns sumPtTracks or sumPtTracks/ptEle according to m_relative
  float getSumPtSquaredTracks();
  
private:

  TLorentzVector m_electron;
  float m_lipEle;
  std::vector<TLorentzVector> m_tracks;
  std::vector<float> m_lips;
  
  float m_extRadius;
  float m_intRadius;
  float m_maxDeltaLip;
  float m_minPtTrack;
  
  bool m_relative;

};

#endif
