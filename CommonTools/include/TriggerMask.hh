//-------------------------------------------------------
// Description:
//    Class to translate the trigger word
//    stored in Conditions tree
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef TRIGGERMASK_H
#define TRIGGERMASK_H

#include "../CommonTools/include/Conditions.hh"
#include "TTree.h"
#include <vector>


class TriggerMask : public Conditions {

public:

  //! constructor
  TriggerMask(TTree *tree=0);
  //! destructor
  virtual ~TriggerMask() {};

  //! get the vector with the requested triggers
  std::vector<int> getBits() { return m_requiredTriggers; }
  //! require a trigger by string
  void requireTrigger(const char* triggerString);
  //! reset the list of required triggers
  void clear() { m_requiredTriggers.clear(); }
  
private:

  //! vector of the requested bits
  std::vector<int> m_requiredTriggers;
  //! the tree holding the conditions variables
  TTree* m_tree;

};

#endif


  
