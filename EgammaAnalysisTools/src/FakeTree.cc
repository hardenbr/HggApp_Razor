#include "EgammaAnalysisTools/include/FakeTree.hh"

FakeTree::FakeTree(const char *filename) {

  myFile = new TFile(filename,"RECREATE");
  myTree = new TTree("T1","fake tree");

  myTree->Branch("met",          &myMet,          "met/F");
  myTree->Branch("Wmt",          &myWmt,          "Wmt/F");
  myTree->Branch("mee",          &myMee,          "mee/F");
  myTree->Branch("jetPt",        &myJetPt,        "jetPt/F");
  myTree->Branch("denomPt",      &myDenomPt,      "denomPt/F");
  myTree->Branch("deltaPhi",     &myDeltaPhi,     "deltaPhi/F");
}

FakeTree::~FakeTree() {
  delete myFile;
}

void FakeTree::store() {
  myTree->Fill();
}

void FakeTree::save() {
  myFile->cd();
  myTree->Write();
  myFile->Close();
}

void FakeTree::fill(float met, float masst, float me, float jpt, float dpt, float dphi){ 
  myMet      = met;
  myWmt      = masst;
  myMee      = me;
  myJetPt    = jpt;
  myDenomPt  = dpt;
  myDeltaPhi = dphi;
}

