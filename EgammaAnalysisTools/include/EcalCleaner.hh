#ifndef EcalCleaner_h
#define EcalCleaner_h

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"

class EcalCleaner {

public:

  //! constructor
  EcalCleaner();

  //! destructor
  virtual ~EcalCleaner();   

  //! configure from files
  void Configure(const char *configDir);

  //! set event by event observables
  void SetE1(float e1) { m_e1 = e1; }
  void SetE4SwissCross(float e4) { m_e4SwissCross = e4; }
  void SetFiducialFlag(int flag) { m_fiducialFlag = flag; }
  void SetSeedFlag(int word) { m_seedFlag = word; }
  void SetSeedTime(float time) { m_seedTime = time; }
  void SetSeedChi2(float chi2) { m_seedChi2 = chi2; }

  //! get output of the selector
  bool output();

  //! display the electron efficiency
  void displayEfficiencies();

private:

  float m_e1;
  float m_e4SwissCross;
  int m_fiducialFlag;
  int m_seedFlag;
  float m_seedTime;
  float m_seedChi2;

  //! contains the selection cuts
  Selection* m_selection;

  //! counters for the efficiencies display, based on cluster candidates
  Counters m_electronCounter;

};

#endif
