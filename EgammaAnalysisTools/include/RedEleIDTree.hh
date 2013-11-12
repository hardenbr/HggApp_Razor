#ifndef RedEleIDTree_H
#define RedEleIDTree_H

#include <TFile.h>
#include <TTree.h>

class RedEleIDTree {
public:
  
  RedEleIDTree(const char *filename);
  ~RedEleIDTree();

  //! add the electron attributes (see below)
  void addAttributesSignal();
  void addAttributesBackground();
  //! add the splitting categories (see below)
  void addCategories();
  void addMore();
  void addTrackMomenta();
  void addElectronIdBits();
  void addDenominatorFakeBits();
  void addIsolations();
  void addGamma();
  //! add run,lumi, event number (for data)
  void addRunInfos();

  //! fill the tree with electron id variables
  //! minimal set
  void fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float See, float Spp, float fbrem, int nbrems, float pt, float eta, int charge);
  //! needed for HZZ BDT
  void fillVariables(float eleEoPout, float EoPout, float EoP, float HoE, float Deta, float Dphi, float s9s25, float s1s9, float See, float Spp, float fbrem, 
                     int nbrems, int nHits, float dcot, float dist, float pt, float eta, int charge, float phiwidth, float etawidth,
                     float IoEmIoP, float eledeta, float d0, float ip3d, float ip3ds, int kfhits, int kflayers, float kfchi2, float e1x5e5x5, int ecaldriven, bool matchConv,
                     bool iseb, bool isee);
  //! additional needed for HWW BDT
  void fillVariables2(float detacalo, float dphicalo, float sep, float dz, float gsfchi2, float emaxovere, float etopovere, float ebottomovere, float eleftovere, float erightovere,
                      float e2ndovere, float e2x5rightovere, float e2x5leftovere, float e2x5topevere, float e2x5bottomovere, 
                      float e2x5maxovere, float e1x5overe, float e2x2overe, float e3x3overe, float e5x5overe, float r9, float nclu,
                      float phi, float scenergy, float scrawenergy, float scesenergy, float eseedopin);

  //! fill ECAL cluster informations
  void fillCluterInfos(float scEt, float scEta, float scPhi, float EtaSeed, float PhiSeed, float ESeed, float IEtaSeed, float IPhiSeed, float EtaCrySeed, float PhiCrySeed, float IEtaCrySeed, float IPhiCrySeed);

  //! fill the tree with isolation variables
  void fillIsolations(float tkIso[2], float ecalIso[2], float hcalIso[2],
                      float combPFiso,
                      float chaPFiso[8], float neuPFiso[8], float phoPFiso[8]);

  //! fill the electron ID bits
  void fillCutBasedIDBits(int CutBasedId[6], int CutBasedIdOnlyID[6], int CutBasedIdOnlyIso[6], int CutBasedIdOnlyConv[6]);
  void fillLHBasedIDBits(int LHBasedId[5], int LHBasedIdOnlyID[5], int LHBasedIdOnlyIso[5], int LHBasedIdOnlyConv[5]);
  void fillLHBasedPFIsoIDBits(int LHBasedPFIsoId[5], int LHBasedPFIsoIdOnlyID[5], int LHBasedPFIsoIdOnlyIso[5], int LHBasedPFIsoIdOnlyConv[5]);
  void fillCiCBasedIDBits(int cic[5]);
  void fillFakeRateDenomBits(float leadJetPt, bool isDenom, bool isDenomTrigger);
  void fillBDTBasedIDBits(int isBDTOnlyId);

  //! fill electron attributes + z mass for the tag and probe
  //! note: when both electrons from Z are probes, the same Z mass is repeated
  void fillAttributesSignal(float zmass, int zeeDec, float genenergy, float genenergystatus1, float genenergystatus3);
  //! fill electron attributes + other quantities for background tag and probe
  void fillAttributesBackground(float dphi, float invmass, float met, float pth);
  //! fill the splitting categories of the PDFs
  void fillCategories(int iecal, int iptbin, int iclass, int nbr);
  void fillMore(float nVtx, float rho, float bdthww[2], float newbdthww[4], float bdthzz[4], float pfmva, float like);
  void fillTrackMomenta(float pcomb, float pmodegsf, float pmeangsf, float pkf, float pterrorgsf, float pterrorkf, float perror);
  void fillGamma(float atg, float aeg, float ahg, int ig);
  //! fill the run,lumi, event number, mc match
  void fillRunInfos(int run, int lumi, int event, int npu[3], int mcmatch);   

  void store();
  void save();

private:

  float myEEleoPout, myEseedoPout, myEoP,myHoE,myDeta,myDphi,mys9s25,mys1s9,mySee,mySpp,myFbrem, myPhiWidth, myEtaWidth, myEseedoPin;
  float myIoEoIoP, myEleDeta, myD0, myIP3d, myIP3dSig, myKFChi2, myE1x5E5x5, myPreShowerOverRaw;
  float myNbrems, myKFHits, myKFLayers, myMissHits;
  bool myMatchConv, myEcalDriven;
  float myDetaCalo, myDphiCalo, mySep, myDZ, myGSFChi2;
  float mySeedEMax,mySeedETop,mySeedEBottom,mySeedELeft,mySeedERight,mySeedE2nd,mySeedE2x5Right,mySeedE2x5Left,mySeedE2x5Top,mySeedE2x5Bottom;
  float mySeedE2x5Max,mySeedE1x5,mySeedE2x2,mySeedE3x3,mySeedE5x5,myR9,myNClusters,myOneMinusSeedE1x5OverE5x5;
  float myDist, myDcot, myCharge;
  float myEta, myPhi, myPt, mySCEnergy, mySCRawEnergy, myEsenergy, myEcalEnergy;
  float myPComb, myPModeGsf, myPMeanGsf, myPKf, myPtErrorGsf, myPtErrorKf, myPError;
  float myNpu[3];
  UInt_t myRun, myLS, myMCMatch;
  ULong64_t myEvent;
  bool myIsEB, myIsEE;
  float myScEt,myScEta,myScPhi,myEtaSeed,myPhiSeed,myESeed,myIEtaSeed,myIPhiSeed,myEtaCrySeed,myPhiCrySeed,myIEtaCrySeed,myIPhiCrySeed;
  bool myEventNumberParity;
  float myZmass, myZDec, myWeight;
  float myGeneratedEnergy,myGeneratedEnergyStatus1,myGeneratedEnergyStatus3;

  int myCutBasedId[6], myCutBasedIdOnlyID[6], myCutBasedIdOnlyIso[6], myCutBasedIdOnlyConv[6];
  int myLHBasedId[5], myLHBasedIdOnlyID[5], myLHBasedIdOnlyIso[5], myLHBasedIdOnlyConv[5];
  int myLHBasedPFIsoId[5], myLHBasedPFIsoIdOnlyID[5], myLHBasedPFIsoIdOnlyIso[5], myLHBasedPFIsoIdOnlyConv[5];
  bool myDenomFake, myPassTriggerDenominator;
  UInt_t myTriggerBit;
  float myLeadJetPt;
  int myCiC[5];
  int myBDTIdOnlyId;

  float myQCDDeltaphi;
  float myQCDInvmass;
  float myQCDMet;
  float myQCDPtHat;

  int myiecal;
  int myiptbin;
  int myiclass;
  int mynbrem;

  float myTrkIso03,myTrkIso04;
  float myEcalIso03,myEcalIso04;
  float myHcalIso03,myHcalIso04;
  float myPFCandCombinedIsoHWW;
  float myPFCandChargedIso[8], myPFCandNeutralIso[8], myPFCandPhotonIso[8];
  float myChargedIso_DR0p0To0p1,myChargedIso_DR0p1To0p2,myChargedIso_DR0p2To0p3,myChargedIso_DR0p3To0p4,myChargedIso_DR0p4To0p5;
  float myGammaIso_DR0p0To0p1,myGammaIso_DR0p1To0p2,myGammaIso_DR0p2To0p3,myGammaIso_DR0p3To0p4,myGammaIso_DR0p4To0p5;
  float myNeutralHadronIso_DR0p0To0p1,myNeutralHadronIso_DR0p1To0p2,myNeutralHadronIso_DR0p2To0p3,myNeutralHadronIso_DR0p3To0p4,myNeutralHadronIso_DR0p4To0p5;

  float myNVtx;
  float myRho;
  float myBdtHww[2], myBdtHzz[4], myNewBdtHww[4], myPFMVA, myLike;

  float myAbsTrackerIsolGammaCand;
  float myAbsEcalIsolGammaCand;
  float myAbsHcalIsolGammaCand;
  int myIsGamma;

  TFile *myFile;
  TTree *myTree;

};

#endif // RedEleIDTree_H
