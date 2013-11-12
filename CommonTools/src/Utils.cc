#include "CommonTools/include/Utils.hh"
#include <math.h>
#include <iostream>

using namespace bits;

std::vector<int> Utils::getTriggers(std::vector<int> requiredTriggers, int firedTrg[4]) {
  std::vector<int> out;

  // unpack the trigger words
  for( int i=0; i<requiredTriggers.size(); i++ ) {
    std::vector<int> in;
    in.push_back(requiredTriggers[i]);
    out.push_back(getTriggersAND(in, firedTrg));
    in.clear();
  }
  return out;
}

bool Utils::getTriggersAND(std::vector<int> requiredTriggers, int firedTrg[4]) {

  // unpack the trigger words
  for( int i=0; i<requiredTriggers.size(); i++ ) {
    
    int block =  requiredTriggers[i]/30;
    int pos = requiredTriggers[i]%30;
    int word = firedTrg[block];

    if ( !( (word >> pos)%2) ) return false;
  }

  return true;

}

bool Utils::getTriggersOR(std::vector<int> requiredTriggers, int firedTrg[4]) {

  // unpack the trigger words
  for( int i=0; i<requiredTriggers.size(); i++ ) {

    int block =  requiredTriggers[i]/30;
    int pos = requiredTriggers[i]%30;
    int word = firedTrg[block];
    
    if ( (word >> pos)%2 ) return true;
  }

  return false;

}

bool Utils:: getL1TriggersOutput(std::vector<int> requiredTriggers, std::vector<int> notRequiredTriggers, int L1FiredTrg[5]) {

  if ( requiredTriggers.size() == 0 && notRequiredTriggers.size() == 0) return true;
  
  bool passRequiredTrigger = false;
  // unpack the trigger words
  for( int i=0; i<requiredTriggers.size(); i++ ) {

    int block =  requiredTriggers[i]/30;
    int pos = requiredTriggers[i]%30;
    int word = L1FiredTrg[block];
    
    if ( (word >> pos)%2 ) {
      passRequiredTrigger = true;
      break;
    }
  }

  if(requiredTriggers.size() > 0 && !passRequiredTrigger) return false;

  bool passNotRequiredTrigger = false;
  // unpack the trigger words
  for( int i=0; i<notRequiredTriggers.size(); i++ ) {

    int block =  notRequiredTriggers[i]/30;
    int pos = notRequiredTriggers[i]%30;
    int word = L1FiredTrg[block];
    
    if ( (word >> pos)%2 ) {
      passNotRequiredTrigger = true;
      break;
    }
  }
  
  if(notRequiredTriggers.size() > 0 && passNotRequiredTrigger) return false;

  return true;
  
}

bool Utils::isInElectronFiducialEta(float eta) {

  return ( fabs(eta) < 1.4442 || // EB
	   (fabs(eta) > 1.560 && fabs(eta) < 2.5 ) // EE
	   );

}

bool Utils::isInECALFiducial(int word) {

  return ( ( (word >> isEE)%2 || (word >> isEB)%2 ) && !((word >> isEBEEGap)%2) );

}

bool Utils::fiducialFlagECAL(int word, ElectronFiducialBit bit) {

  return ( word >> bit )%2;

}

bool Utils::muonIdVal(int word, MuonIdBit bit) {

  return (word >> bit)%2;

}

bool Utils::electronIdVal(int word, ElectronIdBit bit  ) {
  // for each electron ID type, 3 bits are stored:
  // i            j        k
  // notused      iso      ID
  //  return ( ((word >> 3*bit) & 0b001) >> 0 )%2; // commented: not compiling on MacOS X
  return true;
}

bool Utils::isolVal(int word, ElectronIdBit bit  ) {
  // for each electron ID type, 3 bits are stored:
  // i            j        k
  // notused      iso      ID
  //  return ( ((word >> 3*bit) & 0b010) >> 1 )%2; // commented: not compiling on MacOS X
  return true;
}

bool Utils::electronRecoType(int word, ElectronRecoBit bit) {

  return (word >> bit)%2;

}

bool Utils::electronEnergyCorrectionType(int word, ElectronEnergyCorrectionBit bit) {

  return (word >> bit)%2;
  
}

bool Utils::jetIdVal(int word, bits::Version version, bits:: Quality quality) {

  int bit = version*3+quality;
  return (word >> bit)%2;

}
