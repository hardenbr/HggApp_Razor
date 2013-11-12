#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/JetIdBits.h"

class Utils {

public:

  Utils() { }
  //! returns the values of the requested HLT triggers
  std::vector<int> getTriggers(std::vector<int> requiredTriggers, int firedTrg[4]);
  //! returns the AND of the requested HLT triggers
  bool getTriggersAND(std::vector<int> requiredTriggers, int firedTrg[4]);
  //! returns the OR of the requested HLT triggers
  bool getTriggersOR(std::vector<int> requiredTriggers, int firedTrg[4]);
  //! returns the OR of the requested L1 triggers AND NOT the OR of the not requested triggers
  bool getL1TriggersOutput(std::vector<int> requiredTriggers, std::vector<int> notRequiredTriggers, int L1FiredTrg[5]);
  //! returns true if eta belongs to the electron fiducial region:
  //! remove the gap between EB - EE plus some crystal in the bounds (simple, not precise)
  bool isInElectronFiducialEta(float eta);
  //! remove all the ECAL gaps (from electron flags: precise)
  bool isInECALFiducial(int word);
  //! get value of the bit corresponding to ECAL fiducial flags 
  bool fiducialFlagECAL(int word, bits::ElectronFiducialBit bit);
  //! get value of the bit corresponding to certain muon ID 
  bool muonIdVal(int word, bits::MuonIdBit bit);
  //! get value of the eleID bit corresponding to certain electron ID
  bool electronIdVal(int word, bits::ElectronIdBit bit);
  //! get value of the isolation bit corresponding to certain electron ID
  bool isolVal(int word, bits::ElectronIdBit bit);
  //! get value of the bit corresponding to electron reconstruction type
  bool electronRecoType(int word, bits::ElectronRecoBit bit);
  //! get value of the bit corresponding to energy correction type
  bool electronEnergyCorrectionType(int word, bits::ElectronEnergyCorrectionBit bit);
  //! get value of the JetID bit corresponding to certain jet ID Version and Quality
  bool jetIdVal(int word, bits::Version version, bits:: Quality quality);

protected:

};

#endif

