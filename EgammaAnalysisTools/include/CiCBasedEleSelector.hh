#ifndef CiCBasedEleSelector_h
#define CiCBasedEleSelector_h

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"
#include "EgammaAnalysisTools/include/EcalCleaner.hh"

class CiCBasedEleSelector {

public:

  //! constructor
  CiCBasedEleSelector();
  //! destructor
  virtual ~CiCBasedEleSelector();   
  //! configure from files (class dependent)
  void Configure(std::string type, bool useEtBins, bool specialCategories, int version);
  //! configure the ECAL cleaner
  void ConfigureEcalCleaner(const char *configDir);
  //! get the non-class dependent selection (EB, EE)
  Selection* GetSelectionNoClass(const char *EcalSubDet);
  //! set event by event observables
  void SetSCEt(float SCEt) { m_SCEt = SCEt;  }
  void SetSCEta(float SCEta) { m_SCEta = SCEta;  }
  void SetHOverE(float HOverE) { m_HOverE = HOverE; }
  void SetDEta(float DEta) { m_DEta = DEta;  }
  void SetDPhiIn(float DPhiIn) { m_DPhiIn = DPhiIn;  }
  void SetBremFraction(float BremFraction) { m_BremFraction = BremFraction;  }
  void SetSigmaEtaEta(float SigmaEtaEta) { m_SigmaEtaEta = SigmaEtaEta;  }
  void SetEOverPin(float EOverPin) { m_EOverPin = EOverPin;  }
  void SetESeedOverPin(float ESeedOverPin) { m_ESeedOverPin = ESeedOverPin;  }
  void SetEcalFiducialRegion(int word) { m_fiducialflag = word; }
  void SetRecoFlag(int word) { m_recoflag = word; }
  void SetEcalIsolation(float ecalIso) { m_ecalIso = ecalIso;  }
  void SetTrkIsolation(float trackerIso) { m_trackerIso = trackerIso;  }
  void SetHcalIsolation(float hcalIso) { m_hcalIso = hcalIso;  }
  void SetMissingHits(int missingHits) { m_missingHits = missingHits;  }
  void SetConvDist(float dist) { m_distConv = dist;  }
  void SetConvDcot(float dcot) { m_dcotConv = dcot;  }

  void ElectronClassification() ;

  //! get output of the selector (class dependent)
  bool output();
  //! get output of the selector (class dependent) only EleId
  bool outputEleId(); 
  //! get output of the selector (class dependent) only Iso
  bool outputIso();
  //! get output of the selector (class dependent) only Conv
  bool outputConv();
  //! reset
  bool reset();

  void displayEfficiencies();

  bool compute_eid_cut(float x, float et, float cut_min, float cut_max, bool gtn);

  //! ECAL cleaner (public to set the variables through CiCBasedEleSelector)
  EcalCleaner *m_cleaner;

private:

  bool m_doEcalCleaning;

  float m_SCEt;
  float m_SCEta;
  float m_ESeedOverPin;
  float m_HOverE;
  float m_DEta;
  float m_DPhiIn;
  float m_BremFraction;
  float m_SigmaEtaEta;
  float m_EOverPin;
  float m_ecalIso;
  float m_trackerIso;
  float m_hcalIso;
  float m_distConv;
  float m_dcotConv;
  int m_missingHits;
  int m_fiducialflag;
  int m_recoflag;
  int m_cat;
  int m_ptBin;
  bool m_useClass;
  bool m_classDep;
  bool m_useEtBins;
  bool m_specialCategories;

  int m_version;
  int m_eIDLevel;

  //! counters for the efficiencies display, based on electron candidates
  Counters m_electronCounter;

};

#endif
