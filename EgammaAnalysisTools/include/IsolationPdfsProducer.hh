//---------------------------------------
// Description:
// Create the isolation distributions from Z->ee Tag and probe
//---------------------------------------

#ifndef LHPDFSPRODUCER_H
#define LHPDFSPRODUCER_H

#include <vector>

#include "TH1F.h"
#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/EgammaBase.h"

class IsolationPdfsProducer : public EgammaBase {

public:
  //! constructor
  IsolationPdfsProducer(TTree *tree=0);
  //! destructor
  virtual ~IsolationPdfsProducer();
  //! loop over events doing signal pdfs with Zee 
  void LoopZTagAndProbe(const char *treefilesuffix);
  //! loop over events doing bkg pdfs on QCD
  void LoopQCD();
  //! loop over events doing bkg pdfs on W+jets
  void LoopWjets();
  //! set the list of the required triggers
  void requireTrigger(std::vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
  //! save the pdfs in a ROOT file
  void saveHistos(const char *filename);

private:
  
  //! book the histograms of the PDFs
  void bookHistos();

  //! contains the class-dependent electron ID cuts
  std::vector<Selection*> m_EgammaCutIDSelection;

  //! the tag and the probe
  int electrons[2];
  int muons[2];

  //! the required trigger bits
  std::vector<int> m_requiredTriggers;

  //! the configurable selection
  Selection *m_selection;

  // ---------- monitoring histograms ------------
  TH1F *m_Zmass;

  /// Electrons: not splitted
  /// histo[ecalsubdet]
  TH1F *tkSumPtRel[2];
  TH1F *ecalSumEtRel[2];
  TH1F *hcalSumEtRel[2];

  TH1F *ecalRecHitsSumEtRel[2];

};

#endif
