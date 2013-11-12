///---------------------------------------
// Description:
// Class to estimate the performances of Likelihood Electron ID
//---------------------------------------

#ifndef LIKELIHOODANALYSIS_H
#define LIKELIHOODANALYSIS_H

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/CiCBasedEleSelector.hh"
#include "EgammaAnalysisTools/include/Egamma.h"
#include <TLorentzVector.h>

#include <TTree.h>
#include "include/FakeTree.hh"

class LikelihoodAnalysis : public Egamma {

public:
  //! constructor
  LikelihoodAnalysis(TTree *tree=0);
  //! destructor
  virtual ~LikelihoodAnalysis();
  //! find the cut on LH that equals the electron ID with cuts (standard egamma)
  float findEquivalentLHCut();
  //! find the cut on LH that gives a wanted efficiency
  float findEquivalentLHCut(float wantEfficiency);
  //! reproduce the egamma cut based ID (for debugging)
  void reproduceEgammaCutID();
  //! produce the ID efficiency eta/pT distributions
  void estimateIDEfficiency(const char *outname="job0");
  //! produce the mis-ID eta/pT distributions
  void estimateFakeRate(const char *outname="job0");
  //! produce the mis-ID eta/pT distributions from QCD di-jets
  void estimateFakeRateQCD(const char *outname="job0");
  //! produce the mis-ID eta/pT distributions for WW studies
  void estimateFakeRateForHToWW_QCD(const char *outname="job0");
  void estimateFakeRateForHToWW_EGAMMA(const char *outname="job0");

private:

  //! apply the custom offline electron ID
  void isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput, bool applyBDTIdNotCutbased);

  //! apply the custom offline electron ID
  void isEleID(CiCBasedEleSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);  

  float SigmaiEiE(int electron);
  float SigmaiPiP(int electron);

  bool isDenomFake_HwwEgamma(int theEle);
  bool isDenomFake_smurfs(int theEle);

  CutBasedEleIDSelector EgammaCutBasedIDHWW;

  std::vector<CutBasedEleIDSelector> EgammaCutBasedID;
  std::vector<CiCBasedEleSelector>   EgammaCiCBasedID;
  std::vector<CutBasedEleIDSelector> EgammaLHBasedID;
  std::vector<CutBasedEleIDSelector> EgammaBdtBasedID;

  std::vector<std::string> EgammaCutBasedIDWPs;
  std::vector<std::string> EgammaCiCBasedIDWPs;
  std::vector<std::string> EgammaLHBasedIDWPs;
  std::vector<std::string> EgammaBdtBasedIDWPs;

  ElectronLikelihood *LH;
  ElectronIDMVA *fMVA;

  //! contains the class-dependent electron ID cuts
  std::vector<Selection*> m_EgammaCutIDSelection;
  
  bool _isData;

  // counters
  Counters myCounter;

  //! reduced trees      
  FakeTree *myOutTree;
  FakeTree *myOutTreePassed;
};

#endif
