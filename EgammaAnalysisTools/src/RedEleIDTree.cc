#include "EgammaAnalysisTools/include/RedEleIDTree.hh"

RedEleIDTree::RedEleIDTree(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myFile->mkdir("eleIDdir");
  myTree = new TTree("T1","eleID tree");

  // electron kinematics
  myTree->Branch("pt",              &myPt,              "pt/F");
  myTree->Branch("eta",             &myEta,             "eta/F");
  myTree->Branch("phi",             &myPhi,             "phi/F");
  myTree->Branch("Charge",          &myCharge,          "Charge/F");
  myTree->Branch("IsEB",            &myIsEB,            "IsEB/O");
  myTree->Branch("IsEE",            &myIsEE,            "IsEE/O");

  // ele ID vars
  myTree->Branch("EEleoPout",       &myEEleoPout,       "EEleoPout/F");
  myTree->Branch("ESeedoPout",      &myEseedoPout,      "ESeedoPout/F");
  myTree->Branch("EoP",             &myEoP,             "EoP/F");
  myTree->Branch("EoPin",           &myEseedoPin,       "EoPin/F");
  myTree->Branch("IoEmIoP",         &myIoEoIoP,         "IoEmIoP/F");
  myTree->Branch("HoE",             &myHoE,             "HoE/F");
  myTree->Branch("eledeta",         &myEleDeta,         "eledeta/F");
  myTree->Branch("deta",            &myDeta,            "deta/F");
  myTree->Branch("dphi",            &myDphi,            "dphi/F");
  myTree->Branch("detacalo",        &myDetaCalo,        "detacalo/F");
  myTree->Branch("dphicalo",        &myDphiCalo,        "dphicalo/F");
  myTree->Branch("fbrem",           &myFbrem,           "fbrem/F");
  myTree->Branch("nbrems",          &myNbrems,          "nbrems/F");
  myTree->Branch("missHits",        &myMissHits,        "missHits/F");
  myTree->Branch("dist",            &myDist,            "dist/F");
  myTree->Branch("dcot",            &myDcot,            "dcot/F");
  myTree->Branch("d0",              &myD0,              "d0/F");
  myTree->Branch("dz",              &myDZ,              "dz/F");
  myTree->Branch("ip3d",            &myIP3d,            "ip3d/F");
  myTree->Branch("ip3ds",           &myIP3dSig,         "ip3ds/F");
  myTree->Branch("kfhits",          &myKFHits,          "kfhits/F");
  myTree->Branch("kflayers",        &myKFLayers,        "kflayers/F");
  myTree->Branch("kfchi2",          &myKFChi2,          "kfchi2/F");
  myTree->Branch("gsfchi2",         &myGSFChi2,         "gsfchi2/F");

  // cluster shapes
  myTree->Branch("s9s25",           &mys9s25,           "s9s25/F");
  myTree->Branch("phiwidth",        &myPhiWidth,        "phiwidth/F");
  myTree->Branch("etawidth",        &myEtaWidth,        "etawidth/F");
  myTree->Branch("see",             &mySee,             "see/F");
  myTree->Branch("sep",             &mySep,             "sep/F");
  myTree->Branch("spp",             &mySpp,             "spp/F");
  myTree->Branch("NClusters",       &myNClusters,       "NClusters/F");
  myTree->Branch("e1x5e5x5",        &myE1x5E5x5,        "e1x5e5x5/F");
  // cluster positions/energies
  myTree->Branch("scEt",            &myScEt,            "scEt/F");
  myTree->Branch("scEta",           &myScEta,           "scEta/F");
  myTree->Branch("scPhi",           &myScPhi,           "scPhi/F");
  myTree->Branch("SCRawEnergy",     &mySCRawEnergy,     "SCRawEnergy/F"); 
  myTree->Branch("ecalenergy",      &myEcalEnergy,      "ecalenergy/F"); // this is dummy: not in vecbos
  myTree->Branch("esenergy",        &myEsenergy,        "esenergy/F");
  myTree->Branch("PreShowerOverRaw", &myPreShowerOverRaw, "PreShowerOverRaw/F");
  // seed basic cluster
  myTree->Branch("EtaSeed",         &myEtaSeed,         "EtaSeed/F");
  myTree->Branch("PhiSeed",         &myPhiSeed,         "PhiSeed/F");
  myTree->Branch("ESeed",           &myESeed,           "ESeed/F");
  myTree->Branch("IEtaSeed",         &myIEtaSeed,       "IEtaSeed/F");
  myTree->Branch("IPhiSeed",         &myIPhiSeed,       "IPhiSeed/F");
  // seed crystal
  myTree->Branch("EtaCrySeed",       &myEtaCrySeed,     "EtaCrySeed/F");
  myTree->Branch("PhiCrySeed",       &myPhiCrySeed,     "PhiCrySeed/F");
  myTree->Branch("IEtaCrySeed",      &myIEtaCrySeed,    "IEtaCrySeed/F");
  myTree->Branch("IPhiCrySeed",      &myIPhiCrySeed,    "IPhiCrySeed/F");

  // cluster shapes
  myTree->Branch("EMaxSeed",   &mySeedEMax,   "EMaxSeed/F"); 
  myTree->Branch("ETopSeed",   &mySeedETop,   "ETopSeed/F"); 
  myTree->Branch("EBottomSeed",&mySeedEBottom,"EBottomSeed/F"); 
  myTree->Branch("ELeftSeed",  &mySeedELeft,  "ELeftSeed/F"); 
  myTree->Branch("ERightSeed", &mySeedERight, "ERightSeed/F"); 
  myTree->Branch("E2ndSeed",   &mySeedE2nd,   "E2ndSeed/F"); 
  myTree->Branch("E2x5RightSeed",  &mySeedE2x5Right, "E2x5RightSeed/F"); 
  myTree->Branch("E2x5LeftSeed",   &mySeedE2x5Left,  "E2x5LeftSeed/F"); 
  myTree->Branch("E2x5TopSeed",&mySeedE2x5Top,"E2x5TopSeed/F"); 
  myTree->Branch("E2x5BottomSeed", &mySeedE2x5Bottom, "E2x5BottomSeed/F"); 
  myTree->Branch("E2x5MaxSeed",&mySeedE2x5Max,"E2x5MaxSeed/F"); 
  myTree->Branch("E1x5Seed",   &mySeedE1x5,   "E1x5Seed/F"); 
  myTree->Branch("E2x2Seed",   &mySeedE2x2,   "E2x2Seed/F"); 
  myTree->Branch("E3x3Seed",   &mySeedE3x3,   "E3x3Seed/F"); 
  myTree->Branch("E5x5Seed",   &mySeedE5x5,   "E5x5Seed/F"); 
  myTree->Branch("OneMinusSeedE1x5OverE5x5", &myOneMinusSeedE1x5OverE5x5, "OneMinusSeedE1x5OverE5x5/F");
  myTree->Branch("R9",              &myR9,              "R9/F");
  myTree->Branch("matchConv",       &myMatchConv,       "matchConv/O");
  myTree->Branch("ecaldriven",      &myEcalDriven,      "ecaldriven/O");
}

RedEleIDTree::~RedEleIDTree() {
  delete myFile;
}

void RedEleIDTree::addAttributesSignal() {
  myTree->Branch("mass",       &myZmass,       "mass/F");
  myTree->Branch("zdec",       &myZDec,        "zdec/F");
  myTree->Branch("GeneratedEnergy",         &myGeneratedEnergy,         "GeneratedEnergy/F");
  myTree->Branch("GeneratedEnergyStatus1",  &myGeneratedEnergyStatus1,  "GeneratedEnergyStatus1/F");
  myTree->Branch("GeneratedEnergyStatus3",  &myGeneratedEnergyStatus3,  "GeneratedEnergyStatus3/F");
}

void RedEleIDTree::addElectronIdBits() {

  myTree->Branch("WP95",         &myCutBasedId[0],         "WP95/I");
  myTree->Branch("WP90",         &myCutBasedId[1],         "WP90/I");
  myTree->Branch("WP85",         &myCutBasedId[2],         "WP85/I");
  myTree->Branch("WP80",         &myCutBasedId[3],         "WP80/I");
  myTree->Branch("WP70",         &myCutBasedId[4],         "WP70/I");
  myTree->Branch("WPSmurf",      &myCutBasedId[5],         "WPSmurf/I");
  myTree->Branch("CutBasedIdOlyID",    myCutBasedIdOnlyID,   "CutBasedIdOnlyID[6]/I");
  myTree->Branch("CutBasedIdOnlyIso",  myCutBasedIdOnlyIso,  "CutBasedIdOnlyIso[6]/I");
  myTree->Branch("CutBasedIdOnlyConv", myCutBasedIdOnlyConv, "CutBasedIdOnlyConv[6]/I");
  myTree->Branch("cic", myCiC, "cic[5]/I");
  myTree->Branch("BDTIdOnlyId",        &myBDTIdOnlyId,           "BDTIdOnlyId/I");
}

void RedEleIDTree::addDenominatorFakeBits() {
  myTree->Branch("DenomFake",          &myDenomFake,             "DenomFake/O");
  myTree->Branch("PassTriggerDenominator", &myPassTriggerDenominator, "PassTriggerDenominator/O");
  myTree->Branch("leadJetPt",          &myLeadJetPt,             "leadJetPt/F");
  myTree->Branch("triggerBit",         &myTriggerBit,            "TriggerBit/i"); // dummy
}

void RedEleIDTree::addRunInfos() {
  myTree->Branch("run",     &myRun,     "run/i");
  myTree->Branch("lumi",    &myLS,      "lumi/i");
  myTree->Branch("event",   &myEvent,   "event/l");
  myTree->Branch("EventNumberParity",  &myEventNumberParity,  "EventNumberParity/O");
  myTree->Branch("weight", &myWeight, "weight/F"); // dummy=1
  myTree->Branch("npu",      myNpu,     "npu[3]/F");
  myTree->Branch("mcmatch", &myMCMatch, "mcmatch/O");
}

void RedEleIDTree::addAttributesBackground() {

  myTree->Branch("qcdDeltaphi",    &myQCDDeltaphi,    "qcdDeltaphi/F");
  myTree->Branch("qcdInvmass",     &myQCDInvmass,     "qcdInvmass/F");
  myTree->Branch("qcdMet",         &myQCDMet,         "qcdMet/F");
  myTree->Branch("qcdPtHat",       &myQCDPtHat,       "qcdPtHat/F");
}

void RedEleIDTree::addCategories() {

  myTree->Branch("iecal",      &myiecal,      "iecal/I");
  myTree->Branch("iptbin",     &myiptbin,     "iptbin/I");
  myTree->Branch("iclass",     &myiclass,     "iclass/I");
  myTree->Branch("nbrem",      &mynbrem,      "nbrem/I");
}

void RedEleIDTree::addIsolations() {
  myTree->Branch("trkIso03",  &myTrkIso03,    "trkIso03/F");
  myTree->Branch("ecalIso03", &myEcalIso03,   "ecalIso03/F");
  myTree->Branch("hcalIso03", &myHcalIso03,   "hcalIso03/F");
  myTree->Branch("trkIso04",  &myTrkIso04,    "trkIso04/F");
  myTree->Branch("ecalIso04", &myEcalIso04,   "ecalIso04/F");
  myTree->Branch("hcalIso04", &myHcalIso04,   "hcalIso04/F");
  myTree->Branch("combPFIsoHWW", &myPFCandCombinedIsoHWW, "combPFIsoHWW/F");
  myTree->Branch("chaPFIso",     myPFCandChargedIso,     "chPFIso[8]/F");
  myTree->Branch("neuPFIso",     myPFCandNeutralIso,     "neuPFIso[8]/F");
  myTree->Branch("phoPFIso",     myPFCandPhotonIso,      "phoPFIso[8]/F");
  // 0.1 rings (computed from above full cone isolations)
  myTree->Branch("ChargedIso_DR0p0To0p1",  &myChargedIso_DR0p0To0p1,   "ChargedIso_DR0p0To0p1/F");
  myTree->Branch("ChargedIso_DR0p1To0p2",  &myChargedIso_DR0p1To0p2,   "ChargedIso_DR1p0To0p2/F");
  myTree->Branch("ChargedIso_DR0p2To0p3",  &myChargedIso_DR0p2To0p3,   "ChargedIso_DR2p0To0p3/F");
  myTree->Branch("ChargedIso_DR0p3To0p4",  &myChargedIso_DR0p3To0p4,   "ChargedIso_DR3p0To0p4/F");
  myTree->Branch("ChargedIso_DR0p4To0p5",  &myChargedIso_DR0p4To0p5,   "ChargedIso_DR4p0To0p5/F");
  myTree->Branch("GammaIso_DR0p0To0p1",  &myGammaIso_DR0p0To0p1,   "GammaIso_DR0p0To0p1/F");
  myTree->Branch("GammaIso_DR0p1To0p2",  &myGammaIso_DR0p1To0p2,   "GammaIso_DR1p0To0p2/F");
  myTree->Branch("GammaIso_DR0p2To0p3",  &myGammaIso_DR0p2To0p3,   "GammaIso_DR2p0To0p3/F");
  myTree->Branch("GammaIso_DR0p3To0p4",  &myGammaIso_DR0p3To0p4,   "GammaIso_DR3p0To0p4/F");
  myTree->Branch("GammaIso_DR0p4To0p5",  &myGammaIso_DR0p4To0p5,   "GammaIso_DR4p0To0p5/F");
  myTree->Branch("NeutralHadronIso_DR0p0To0p1",  &myNeutralHadronIso_DR0p0To0p1,   "NeutralHadronIso_DR0p0To0p1/F");
  myTree->Branch("NeutralHadronIso_DR0p1To0p2",  &myNeutralHadronIso_DR0p1To0p2,   "NeutralHadronIso_DR1p0To0p2/F");
  myTree->Branch("NeutralHadronIso_DR0p2To0p3",  &myNeutralHadronIso_DR0p2To0p3,   "NeutralHadronIso_DR2p0To0p3/F");
  myTree->Branch("NeutralHadronIso_DR0p3To0p4",  &myNeutralHadronIso_DR0p3To0p4,   "NeutralHadronIso_DR3p0To0p4/F");
  myTree->Branch("NeutralHadronIso_DR0p4To0p5",  &myNeutralHadronIso_DR0p4To0p5,   "NeutralHadronIso_DR4p0To0p5/F");
}

void RedEleIDTree::addMore() {
  myTree->Branch("bdthww",     myBdtHww,    "bdthww[2]/F");
  myTree->Branch("newbdthww",  myNewBdtHww, "newbdthww[4]/F");
  myTree->Branch("bdthzz",     myBdtHzz,    "bdthzz[4]/F");
  myTree->Branch("lh",       &myLike,   "lh/F");
  myTree->Branch("PFMVA",    &myPFMVA,  "PFMVA/F");
  myTree->Branch("vertices", &myNVtx, "vertices/F");
  myTree->Branch("rho",      &myRho,  "rho/F");
}

void RedEleIDTree::addTrackMomenta() {
  myTree->Branch("pcomb",    &myPComb,    "pcomb/F");
  myTree->Branch("pmodegsf", &myPModeGsf, "pmodegsf/F");
  myTree->Branch("pmeangsf", &myPMeanGsf, "pmeangsf/F");
  myTree->Branch("pmeankf",  &myPKf,      "pmeankf/F");
  myTree->Branch("pterrorgsf", &myPtErrorGsf, "pterrorgsf/F");
  myTree->Branch("pterrorkf",  &myPtErrorKf,  "pterrorkf/F");
  myTree->Branch("perror",  &myPError,  "perror/F");
}

void RedEleIDTree::addGamma() {

  myTree->Branch("absTrackerIsolGammaCand",&myAbsTrackerIsolGammaCand,"absTrackerIsolGammaCand/F");
  myTree->Branch("absEcalIsolGammaCand",   &myAbsEcalIsolGammaCand,   "absEcalIsolGammaCand/F");
  myTree->Branch("absHcalIsolGammaCand",   &myAbsHcalIsolGammaCand,   "absHcalIsolGammaCand/F");
  myTree->Branch("isGamma",                &myIsGamma,                "isGamma/I");
}

void RedEleIDTree::store() {

  myTree->Fill();
}

void RedEleIDTree::save() {

  myFile->cd("eleIDdir");
  myTree->Write();
  myFile->Close();
}

void RedEleIDTree::fillVariables(float EoPout, float EoP, float HoE, float DEta, float DPhi, float s9s25, float See, float Spp, float fbrem, int nbrems, float pt, float eta, int charge) {
  myEseedoPout=EoPout;
  myEoP=EoP;
  myHoE=HoE;
  myDeta=DEta;
  myDphi=DPhi;
  mys9s25=s9s25;
  mySee=See;
  mySpp=Spp;
  myFbrem=fbrem;
  myNbrems=float(nbrems);
  myPt=pt;
  myEta=eta;
  myCharge=float(charge);
}

void RedEleIDTree::fillVariables(float eleEoPout, float EseedoPout, float EoP, float HoE, float Deta, float Dphi, float s9s25, float s1s9, float See, float Spp, float fbrem, 
                                 int nbrems, int nHits, float dcot, float dist, float pt, float eta, int charge, float phiwidth, float etawidth,
                                 float IoEmIoP, float eledeta, float d0, float ip3d, float ip3ds, int kfhits, int kflayers, float kfchi2, float e1x5e5x5, int ecaldriven, bool matchConv, 
                                 bool iseb, bool isee) {
  myEEleoPout=eleEoPout;
  myEseedoPout=EseedoPout;
  myEoP=EoP;
  myHoE=HoE;
  myDeta=Deta;
  myDphi=Dphi;
  mys9s25=s9s25;
  mys1s9=s1s9;
  mySee=See;
  mySpp=Spp;
  myFbrem=fbrem;
  myNbrems=float(nbrems);
  myMissHits=float(nHits);
  myDist=dist;
  myDcot=dcot;
  myPt=pt;
  myEta=eta;
  myCharge=charge;
  myPhiWidth=phiwidth;
  myEtaWidth=etawidth;
  myIoEoIoP=IoEmIoP;
  myEleDeta=eledeta;
  myD0=d0;
  myIP3d=ip3d;
  myIP3dSig=ip3ds;
  myKFHits=float(kfhits);
  myKFLayers=float(kflayers);
  myKFChi2=kfchi2;
  myE1x5E5x5=e1x5e5x5;
  myEcalDriven=(ecaldriven==1) ? true : false;
  myMatchConv=matchConv;
  myIsEB=iseb;
  myIsEE=isee;
}

void RedEleIDTree::fillVariables2(float detacalo, float dphicalo, float sep, float dz, float gsfchi2, float emax, float etop, float ebottom, float eleft, float eright,
                                  float e2nd, float e2x5right, float e2x5left, float e2x5top, float e2x5bottom, 
                                  float e2x5max, float e1x5, float e2x2, float e3x3, float e5x5, float r9, float nclu,
                                  float phi, float scenergy, float scrawenergy, float scesenergy, float eseedopin) {
  myDetaCalo=detacalo;
  myDphiCalo=dphicalo;
  mySep=sep;
  myDZ=dz;
  myGSFChi2=gsfchi2;
  mySeedEMax=emax;
  mySeedETop=etop;
  mySeedEBottom=ebottom;
  mySeedELeft=eleft;
  mySeedERight=eright;
  mySeedE2nd=e2nd;
  mySeedE2x5Right=e2x5right;
  mySeedE2x5Left=e2x5left;
  mySeedE2x5Top=e2x5top;
  mySeedE2x5Bottom=e2x5bottom;
  mySeedE2x5Max=e2x5max;
  mySeedE1x5=e1x5;
  mySeedE2x2=e2x2;
  mySeedE3x3=e3x3;
  mySeedE5x5=e5x5;
  myOneMinusSeedE1x5OverE5x5=1-e1x5/e5x5;
  myR9=r9;
  myNClusters=nclu;
  myPhi=phi;
  mySCEnergy=scenergy;
  mySCRawEnergy=scrawenergy;
  myEcalEnergy=-999.;
  myEsenergy=scesenergy;
  myEseedoPin=eseedopin;
  myPreShowerOverRaw=scesenergy/scrawenergy;
}

void RedEleIDTree::fillCluterInfos(float scEt, float scEta, float scPhi, float EtaSeed, float PhiSeed, float ESeed, float IEtaSeed, float IPhiSeed, 
                                   float EtaCrySeed, float PhiCrySeed, float IEtaCrySeed, float IPhiCrySeed) {
  myScEt=scEt;
  myScEta=scEta;
  myScPhi=scPhi;
  myEtaSeed=EtaSeed;
  myPhiSeed=PhiSeed;
  myESeed=ESeed;
  myIEtaSeed=IEtaSeed;
  myIPhiSeed=IPhiSeed;
  myEtaCrySeed=EtaCrySeed;
  myPhiCrySeed=PhiCrySeed;
  myIEtaCrySeed=IEtaCrySeed;
  myIPhiCrySeed=IPhiCrySeed;
}

void RedEleIDTree::fillIsolations(float trkIso[2], float ecalIso[2], float hcalIso[2],
                                  float combPFiso,
                                  float chaPFiso[8], float neuPFiso[8], float phoPFiso[8]) {
  myTrkIso03=trkIso[0];
  myEcalIso03=ecalIso[0];
  myHcalIso03=hcalIso[0];
  myTrkIso04=trkIso[1];
  myEcalIso04=ecalIso[1];
  myHcalIso04=hcalIso[1];
  myPFCandCombinedIsoHWW=combPFiso;
  for(int i=0;i<8;i++) {
    myPFCandChargedIso[i]=chaPFiso[i];
    myPFCandNeutralIso[i]=neuPFiso[i];
    myPFCandPhotonIso[i]=phoPFiso[i];
  }
  // rings
  myChargedIso_DR0p0To0p1=chaPFiso[0];
  myChargedIso_DR0p1To0p2=chaPFiso[1]-chaPFiso[0];
  myChargedIso_DR0p2To0p3=chaPFiso[2]-chaPFiso[1];
  myChargedIso_DR0p3To0p4=chaPFiso[3]-chaPFiso[2];
  myChargedIso_DR0p4To0p5=chaPFiso[4]-chaPFiso[3];

  myGammaIso_DR0p0To0p1=phoPFiso[0];
  myGammaIso_DR0p1To0p2=phoPFiso[1]-phoPFiso[0];
  myGammaIso_DR0p2To0p3=phoPFiso[2]-phoPFiso[1];
  myGammaIso_DR0p3To0p4=phoPFiso[3]-phoPFiso[2];
  myGammaIso_DR0p4To0p5=phoPFiso[4]-phoPFiso[3];

  myNeutralHadronIso_DR0p0To0p1=neuPFiso[0];
  myNeutralHadronIso_DR0p1To0p2=neuPFiso[1]-neuPFiso[0];
  myNeutralHadronIso_DR0p2To0p3=neuPFiso[2]-neuPFiso[1];
  myNeutralHadronIso_DR0p3To0p4=neuPFiso[3]-neuPFiso[2];
  myNeutralHadronIso_DR0p4To0p5=neuPFiso[4]-neuPFiso[3];

}

void RedEleIDTree::fillAttributesSignal(float zmass, int zdec, float genenergy, float genenergystatus1, float genenergystatus3) {
  myZmass=zmass;
  myZDec=float(zdec);
  myGeneratedEnergy=genenergy;
  myGeneratedEnergyStatus1=genenergystatus1;
  myGeneratedEnergyStatus3=genenergystatus3;
}

void RedEleIDTree::fillAttributesBackground(float deltaphi, float invmass, float met, float pth) {

  myQCDDeltaphi=deltaphi;
  myQCDInvmass=invmass;
  myQCDMet=met;
  myQCDPtHat=pth;
}

void RedEleIDTree::fillCategories(int iecal, int iptbin, int iclass, int nbr) {

  myiecal=iecal;
  myiptbin=iptbin;
  myiclass=iclass;
  mynbrem=nbr;
}

void RedEleIDTree::fillMore(float nVtx, float rho, float bdthww[2], float newbdthww[4], float bdthzz[4], float pfmva, float like) {
  myNVtx=nVtx;
  myRho=rho;
  for(int i=0; i<2; i++) myBdtHww[i]=bdthww[i];
  for(int i=0; i<4; i++) myNewBdtHww[i]=newbdthww[i];
  for(int i=0; i<4; i++) myBdtHzz[i]=bdthzz[i];
  myPFMVA=pfmva;
  myLike=like;
}

void RedEleIDTree::fillTrackMomenta(float pcomb, float pmodegsf, float pmeangsf, float pkf, float pterrorgsf, float pterrorkf, float perror) {
  myPComb=pcomb;
  myPModeGsf=pmodegsf;
  myPMeanGsf=pmeangsf;
  myPKf=pkf;
  myPtErrorGsf=pterrorgsf;
  myPtErrorKf=pterrorkf;
  myPError=perror;
}

void RedEleIDTree::fillGamma(float atg, float aeg, float ahg, int ig) {

  myAbsTrackerIsolGammaCand = atg;
  myAbsEcalIsolGammaCand    = aeg;
  myAbsHcalIsolGammaCand    = ahg;
  myIsGamma=ig;
}

void RedEleIDTree::fillCutBasedIDBits(int CutBasedId[6], int CutBasedIdOnlyID[6], int CutBasedIdOnlyIso[6], int CutBasedIdOnlyConv[6]) {
  for(int i=0; i<6; i++) {
    myCutBasedId[i] = CutBasedId[i];
    myCutBasedIdOnlyID[i] = CutBasedIdOnlyID[i];
    myCutBasedIdOnlyIso[i] = CutBasedIdOnlyIso[i];
    myCutBasedIdOnlyConv[i] = CutBasedIdOnlyConv[i];
  }
}

void RedEleIDTree::fillLHBasedIDBits(int LHBasedId[5], int LHBasedIdOnlyID[5], int LHBasedIdOnlyIso[5], int LHBasedIdOnlyConv[5]) {
  for(int i=0; i<5; i++) {
    myLHBasedId[i] = LHBasedId[i];
    myLHBasedIdOnlyID[i] = LHBasedIdOnlyID[i];
    myLHBasedIdOnlyIso[i] = LHBasedIdOnlyIso[i];
    myLHBasedIdOnlyConv[i] = LHBasedIdOnlyConv[i];
  }
}

void RedEleIDTree::fillLHBasedPFIsoIDBits(int LHBasedPFIsoId[5], int LHBasedPFIsoIdOnlyID[5], int LHBasedPFIsoIdOnlyIso[5], int LHBasedPFIsoIdOnlyConv[5]) {
  for(int i=0; i<5; i++) {
    myLHBasedPFIsoId[i] = LHBasedPFIsoId[i];
    myLHBasedPFIsoIdOnlyID[i] = LHBasedPFIsoIdOnlyID[i];
    myLHBasedPFIsoIdOnlyIso[i] = LHBasedPFIsoIdOnlyIso[i];
    myLHBasedPFIsoIdOnlyConv[i] = LHBasedPFIsoIdOnlyConv[i];
  }
}

void RedEleIDTree::fillCiCBasedIDBits(int cic[5]) {
  for(int i=0;i<5;i++) myCiC[i]=cic[i];
}

void RedEleIDTree::fillFakeRateDenomBits(float leadJetPt, bool isDenom, bool passTriggerDenominator) {
  myDenomFake = isDenom;
  myPassTriggerDenominator = passTriggerDenominator;
  myLeadJetPt = leadJetPt;
  myTriggerBit=-1;
}

void RedEleIDTree::fillBDTBasedIDBits(int isBDTOnlyId) {
  myBDTIdOnlyId = isBDTOnlyId;
}

void RedEleIDTree::fillRunInfos(int run, int lumi, int event, int npu[3], int mcmatch) {
  myRun = run;
  myLS = lumi;
  myEvent = event;
  myEventNumberParity = (event%2==0) ? true : false;
  myWeight=1.0;
  for(int i=0;i<3;i++) myNpu[i]=float(npu[i]);
  myMCMatch = (mcmatch==1) ? true : false;
}
