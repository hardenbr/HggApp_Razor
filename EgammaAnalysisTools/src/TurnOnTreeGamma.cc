#include "EgammaAnalysisTools/include/TurnOnTreeGamma.hh"

TurnOnTreeGamma::TurnOnTreeGamma(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myFile->mkdir("turnoncurvedir");
  myTreeGamma = new TTree("T1","tree for turn on curves");

  // offline photon kinematics
  myTreeGamma->Branch("ptGamma",   &myPtGamma,   "ptGamma/F");
  myTreeGamma->Branch("etaGamma",  &myEtaGamma,  "etaGamma/F");
  myTreeGamma->Branch("phiGamma",  &myPhiGamma,  "phiGamma/F");

  // photon ID bool
  myTreeGamma->Branch("idGamma",  &myIdGamma,  "idGamma/I");
  
  // match photon - HLT candidate
  myTreeGamma->Branch("hltmatchGamma", &myHltmatchGamma, "hltmatchGamma/I");
}

TurnOnTreeGamma::~TurnOnTreeGamma() {

  delete myFile;
}

void TurnOnTreeGamma::addTandPinfo() {

  myTreeGamma->Branch("mass",       &myZmass,       "mass/F");
  myTreeGamma->Branch("ptEle1",     &myPtEle1,      "ptEle1/F");
  myTreeGamma->Branch("ptEle2",     &myPtEle2,      "ptEle2/F");
  myTreeGamma->Branch("etaEle1",    &myEtaEle1,     "etaEle1/F");
  myTreeGamma->Branch("etaEle2",    &myEtaEle2,     "etaEle2/F");
  myTreeGamma->Branch("tagEle1",    &myTagEle1,     "tagEle1/I");
  myTreeGamma->Branch("tagEle2",    &myTagEle2,     "tagEle2/I");  
  myTreeGamma->Branch("isAZ",       &myIsAZ,        "isAZ/I");
  myTreeGamma->Branch("isAZ_tight", &myIsAZ_tight,  "isAZ_tight/I");
}

void TurnOnTreeGamma::addHLTmatchInfo() {

  myTreeGamma->Branch("passHLT20",  &myPassHLT20,  "passHLT20/I");
  myTreeGamma->Branch("passHLT30",  &myPassHLT30,  "passHLT30/I");
  myTreeGamma->Branch("passHLT50",  &myPassHLT50,  "passHLT50/I");
  myTreeGamma->Branch("passHLT75",  &myPassHLT75,  "passHLT75/I");
  myTreeGamma->Branch("passHLT90",  &myPassHLT90,  "passHLT90/I");
}

void TurnOnTreeGamma::addRunInfos() {

  myTreeGamma->Branch("run",     &myRun,     "run/i");
  myTreeGamma->Branch("lumi",    &myLS,      "lumi/i");
  myTreeGamma->Branch("event",   &myEvent,   "event/l");
}

void TurnOnTreeGamma::addMore() {

  myTreeGamma->Branch("vertices", &myNVtx, "vertices/F");
  myTreeGamma->Branch("rho",      &myRho,  "rho/F");
}

void TurnOnTreeGamma::store() {

  myTreeGamma->Fill();
}

void TurnOnTreeGamma::save() {

  myFile->cd("turnoncurvedir");
  myTreeGamma->Write();
  myFile->Close();
}

void TurnOnTreeGamma::fillVariables(float ptg, float etag, float phig, int idg, int hltmatchg) {

  myPtGamma       = ptg;
  myEtaGamma      = etag;
  myPhiGamma      = phig;
  myIdGamma       = idg;
  myHltmatchGamma = hltmatchg;
}

void TurnOnTreeGamma::fillTandPinfo(float mass, float pt1, float pt2, float eta1, float eta2, int tag1, int tag2, int isaz, int isazt) {

  myZmass      = mass; 
  myPtEle1     = pt1;
  myPtEle2     = pt2;
  myEtaEle1    = eta1;
  myEtaEle2    = eta2;
  myTagEle1    = tag1;
  myTagEle2    = tag2;
  myIsAZ       = isaz;
  myIsAZ_tight = isazt;
}

void TurnOnTreeGamma::fillHLTmatchInfo(int hlt20, int hlt30, int hlt50, int hlt75, int hlt90) {

  myPassHLT20 = hlt20;
  myPassHLT30 = hlt30;
  myPassHLT50 = hlt50;
  myPassHLT75 = hlt75;
  myPassHLT90 = hlt90;
}

void TurnOnTreeGamma::fillRunInfos(int run, int lumi, int event) {

  myRun   = run;
  myLS    = lumi;
  myEvent = event;
}

void TurnOnTreeGamma::fillMore(float nVtx, float rho) {

  myNVtx=nVtx;
  myRho=rho;
}

