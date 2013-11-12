#ifndef LEPTON_ID_BITS_H
#define LEPTON_ID_BITS_H

namespace bits {

enum MuonIdBit { TMLastStationOptimizedLowPtTight=0, TMLastStationOptimizedLowPtLoose=1, 
                 TMOneStationTight=2, TMOneStationLoose=3, 
                 TM2DCompatibilityTight=4, TM2DCompatibilityLoose=5,
                 TMLastStationTight=6, TMLastStationLoose=7,
                 GlobalMuonPromptTight=8,
                 AllArbitrated=9, TrackerMuonArbitrated=10,
                 AllTrackerMuons=11, AllStandAloneMuons=12, AllGlobalMuons=13 };

enum ElectronIdBit { eleIdRobustHighEnergy=0, eleIdRobustTight=1, eleIdRobustLoose=2, eleIdTight=3, eleIdLoose=4,
                     eleIdHyperTight1CIC=5, eleIdSuperTightCIC=6, eleIdTightCIC=7, eleIdMediumCIC=8, eleIdLooseCIC=9,
                     eleIdVeryLooseCIC=10};
 
enum ElectronFiducialBit { isEERingGap=0, isEEDeeGap=1, isEEGap=2, isEBPhiGap=3, isEBEtaGap=4, isEBGap=5, isEBEEGap=6, isGap=7, isEE=8, isEB=9 };

enum ElectronRecoBit { isTrackerDriven=0, isEcalDriven=1 };

enum ElectronEnergyCorrectionBit { isMomentumCorrected=0, isEcalEnergyCorrected=1 };

enum ElectronClassification { UNKNOWNCLASS =-1, GOLDEN, BIGBREM, NARROW, SHOWERING, GAP } ;

}

#endif
