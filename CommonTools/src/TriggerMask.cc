#include "CommonTools/include/TriggerMask.hh"
#include <iostream>

using namespace std;

TriggerMask::TriggerMask(TTree *tree)
  : Conditions(tree) {

  m_tree = tree;

}

void TriggerMask::requireTrigger( const char* triggerString ) {
  
  m_tree->GetEntry(0);

  bool found = false;
  for(unsigned int i=0; i<nameHLT->size(); i++) {

    if( !strcmp (triggerString, nameHLT->at(i).c_str() ) ) {
      m_requiredTriggers.push_back( indexHLT->at(i) ) ;
      found = true;
    }

  }

  if ( !found ) cout << "ERROR! Trigger Path = " << triggerString << " not found in the saved paths" << endl;
  
}

