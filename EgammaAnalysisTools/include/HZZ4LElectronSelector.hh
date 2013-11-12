///---------------------------------------
// Description:
// Class to estimate the performances of Likelihood Electron ID
//---------------------------------------

#ifndef HZZ4LELECTRONSELECTOR_H
#define HZZ4LELECTRONSELECTOR_H

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/CiCBasedEleSelector.hh"
#include "EgammaAnalysisTools/include/Egamma.h"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "EgammaAnalysisTools/include/ElectronIDMVA.h"
#include "EgammaAnalysisTools/include/ElectronIDMVAHZZ.h"
#include <TLorentzVector.h>

class HZZ4LElectronSelector : public Egamma {

public:
  //! constructor
  HZZ4LElectronSelector(TTree *tree=0);
  //! destructor
  virtual ~HZZ4LElectronSelector();
  //! do the loop over the events
  void Loop(const char *treefilesuffix);
  
private:

  //! match the HZZ electrons (charge and angle)
  bool mcMatches(int probe);
  
  //! for the prompt rate for fake rate
  int isDenomFake(int eleIndex);
  int isDenomFake_smurfs(int theEle);
  
  float SigmaiEiE(int electron);
  float SigmaiPiP(int electron);
  
  ElectronLikelihood *LH;

  bool isData_;

  //! the configurable selection
  Selection *m_selection;

  // the tag and the probe
  std::vector<int> electrons;

  /// MVA for electron ID. To be created and initialized from the children classes
  ElectronIDMVA *fMVAHWW, *fMVAHWWWithIso;
  ElectronIDMVAHZZ *fMVAHZZDanV0, *fMVAHZZSiV0, *fMVAHZZSiV1, *fMVAHZZSiDanV2;
  ElectronIDMVAHZZ *fMVAHWWDanV0, *fMVAHWWSiV0, *fMVAHWWSiV1, *fMVAHWWSiDanV2;

};

#endif
