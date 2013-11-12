///---------------------------------------
// Description:
// Class to make trigger turn-on curves
//---------------------------------------

#ifndef TESTTURNONCURVEGAMMA_H
#define TESTTURNONCURVEGAMMA_H

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/Egamma.h"
#include <TLorentzVector.h>
#include <TTree.h>

class TestTurnOnCurveGamma : public Egamma {

public:

  //! constructor
  TestTurnOnCurveGamma(TTree *tree=0);
  //! destructor
  virtual ~TestTurnOnCurveGamma();
  // ! do the analysis
  void measureTurnOnGamma(const char *outname);

private:

  std::vector<CutBasedEleIDSelector> EgammaCutBasedIDLowPt;
  std::vector<std::string> EgammaCutBasedIDLowPtWPs;

  //! apply the custom offline electron ID
  void isEleID(CutBasedEleIDSelector *selector, int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);

  //! apply the analysis photonID
  bool isPhotonID(int theGamma);
  int effectiveAreaRegion(float theEta);

  bool _isData;
};

#endif
