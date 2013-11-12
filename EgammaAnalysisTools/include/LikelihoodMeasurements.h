#ifndef LikelihoodMeasurements_h
#define LikelihoodMeasurements_h

struct LikelihoodMeasurements {
  int subdet; // 0=EB, 1=EE
  float pt;
  float deltaPhi,
    deltaEta,
    eSeedClusterOverPout,
    eSuperClusterOverP,
    hadronicOverEm,
    sigmaIEtaIEta,
    sigmaIPhiIPhi,
    fBrem,
    OneOverEMinusOneOverP;
  int nBremClusters;
};

#endif

