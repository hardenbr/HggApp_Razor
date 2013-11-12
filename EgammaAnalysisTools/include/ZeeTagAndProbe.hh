///---------------------------------------
// Description:
// Class to estimate the performances of Likelihood Electron ID
//---------------------------------------

#ifndef ZTAGANDPROBE_H
#define ZTAGANDPROBE_H

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/CiCBasedEleSelector.hh"
#include "EgammaAnalysisTools/include/Egamma.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "EgammaAnalysisTools/include/ElectronIDMVA.h"
#include "EgammaAnalysisTools/include/ElectronIDMVAHZZ.h"
#include <TLorentzVector.h>

class ZeeTagAndProbe : public Egamma {

public:
  //! constructor
  ZeeTagAndProbe(TTree *tree=0);
  //! destructor
  virtual ~ZeeTagAndProbe();
  //! do the loop over the events
  void Loop(const char *treefilesuffix);
  
private:

  //! configurations
  void configSelection(Selection* selection);

  //! apply the custom offline electron ID
  void isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);

  //! apply the custom offline electron ID
  void isEleID(CiCBasedEleSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);

  //! match the Z->ee MC truth: always match the probe, if tag==-1 then the tag is not matched
  bool mcMatches(int probe, int tag=-1);

  //! match the electron with the generator level one coming from the Z
  float genEleEnergy(int ele, int status=-1, float deltar=0.3);

  //! for the prompt rate for fake rate
  int isDenomFake(int eleIndex);
  int isDenomFake_smurfs(int theEle);
  
  float SigmaiEiE(int electron);
  float SigmaiPiP(int electron);
  
  //! class members
  std::vector<CutBasedEleIDSelector> EgammaCutBasedIDLowPt, EgammaCutBasedIDHighPt;
  std::vector<CutBasedEleIDSelector> EgammaLHBasedIDLowPt, EgammaLHBasedIDHighPt;
  std::vector<CutBasedEleIDSelector> EgammaLHBasedPFIsoIDLowPt, EgammaLHBasedPFIsoIDHighPt;

  std::vector<std::string> EgammaCutBasedIDLowPtWPs, EgammaCutBasedIDHighPtWPs;
  std::vector<std::string> EgammaLHBasedIDLowPtWPs, EgammaLHBasedIDHighPtWPs;
  std::vector<std::string> EgammaLHBasedPFIsoIDLowPtWPs, EgammaLHBasedPFIsoIDHighPtWPs;

  ElectronLikelihood *LH;

  bool isData_;

  //! the configurable selection
  Selection *m_selection;

  // the tag and the probe
  int electrons[2];

  /// MVA for electron ID. To be created and initialized from the children classes
  ElectronIDMVA *fMVAHWW, *fMVAHWWWithIso;
  ElectronIDMVAHZZ *fMVAHZZDanV0, *fMVAHZZSiV0, *fMVAHZZSiV1, *fMVAHZZSiDanV2;
  ElectronIDMVAHZZ *fMVAHWWDanV0, *fMVAHWWSiV0, *fMVAHWWSiV1, *fMVAHWWSiDanV2;

};

#endif
