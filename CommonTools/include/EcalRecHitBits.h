#ifndef ECALRECHIT_BITS_H
#define ECALRECHIT_BITS_H

namespace EcalRecHitBits {

  enum Flags { 
    kGood=0,                   // channel ok, the energy and time measurement are reliable
    kPoorReco,                 // the energy is available from the UncalibRecHit, but approximate (bad shape, large chi2)
    kOutOfTime,                // the energy is available from the UncalibRecHit (sync reco), but the event is out of time
    kFaultyHardware,           // The energy is available from the UncalibRecHit, channel is faulty at some hardware level (e.g. noisy)
    kNoisy,                    // the channel is very noisy
    kPoorCalib,                // the energy is available from the UncalibRecHit, but the calibration of the channel is poor
    kSaturated,                // saturated channel (recovery not tried)
    kLeadingEdgeRecovered,     // saturated channel: energy estimated from the leading edge before saturation
    kNeighboursRecovered,      // saturated/isolated dead: energy estimated from neighbours
    kTowerRecovered,           // channel in TT with no data link, info retrieved from Trigger Primitive
    kFake,                     // the signal in the channel is a fake (e.g. a so-called spike)
    kFakeNeighbours,           // the signal in the channel is a fake and it is detected by looking at the neighbours
    kDead,                     // channel is dead and any recovery fails
    kKilled,                   // MC only flag: the channel is killed in the real detector
    kUnknown                   // to easy the interface with functions returning flags
  };

}

#endif
