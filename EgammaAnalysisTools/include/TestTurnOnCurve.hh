///---------------------------------------
// Description:
// Class to make trigger turn-on curves
//---------------------------------------

#ifndef TESTTURNONCURVE_H
#define TESTTURNONCURVE_H

#include <vector>

#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/Egamma.h"
#include <TLorentzVector.h>
#include <TTree.h>

class TestTurnOnCurve : public Egamma {

public:

  //! constructor
  TestTurnOnCurve(TTree *tree=0);
  //! destructor
  virtual ~TestTurnOnCurve();

  void measureTurnOn(const char *outname);

private:

  bool isDenomTurnOn(int theEle); 

  bool _isData;
};

#endif
