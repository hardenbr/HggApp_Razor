//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 15 01:31:13 2010 by ROOT version 5.22/00a
// from TTree dataset/dataset with sWeights
// found on file: ../results/sPlotsTree/sPlots_tree.root
//////////////////////////////////////////////////////////

#ifndef sPlotsPdfsComparison_h
#define sPlotsPdfsComparison_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

#include <iostream>
#include <vector>

#include "EgammaAnalysisTools/include/ElectronLikelihood.h"

class sPlotsPdfsComparison {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types for data
   Double_t        N_sig_sw;
   Double_t        L_N_sig;
   Double_t        N_qcd_sw;
   Double_t        L_N_qcd;
   Double_t        trackerIso;
   Double_t        ecalJIso;
   Double_t        ecalGTIso;
   Double_t        hcalIso;
   Double_t        combinedIso;
   Double_t        classification;
   Double_t        deta;
   Double_t        dphi;
   Double_t        hoe;
   Double_t        see;
   Double_t        spp;
   Double_t        eop;
   Double_t        esc;
   Double_t        pin;
   Double_t        rho;
   Double_t        fbrem;
   Double_t        nbrem;
   Double_t        met;
   Double_t        tcmet;
   Double_t        pfmet;
   Double_t        mt;
   Double_t        tcmt;
   Double_t        pfmt;
   Double_t        pt;
   Double_t        eta;
   Double_t        phi;
   Double_t        charge;
   Double_t        weight;
   Double_t        nPFJets;
   Double_t        nJets;
   Double_t        event;
   Double_t        isIdWP70;
   Double_t        isIdWP80;
   Double_t        isIdWP85;
   Double_t        isIdWP90;
   Double_t        isIdWP95;
   Double_t        isIsoWP70;
   Double_t        isIsoWP80;
   Double_t        isIsoWP85;
   Double_t        isIsoWP90;
   Double_t        isIsoWP95;
   Double_t        isConvRejWP70;
   Double_t        isConvRejWP80;
   Double_t        isConvRejWP85;
   Double_t        isConvRejWP90;
   Double_t        isConvRejWP95;
   Double_t        isWP70;
   Double_t        isWP80;
   Double_t        isWP85;
   Double_t        isWP90;
   Double_t        isWP95;

   Int_t         f_nPU;
   Float_t       f_trackerIso;
   Float_t       f_ecalJIso;
   Float_t       f_ecalGTIso;
   Float_t       f_hcalIso;
   Float_t       f_combinedIso;
   Float_t       f_classification;
   Float_t       f_deta;
   Float_t       f_dphi;
   Float_t       f_hoe;
   Float_t       f_see;
   Float_t       f_spp;
   Float_t       f_eop;
   Float_t       f_esc;
   Float_t       f_pin;
   Float_t       f_rho;
   Float_t       f_fbrem;
   Int_t         f_nbrem;
   Float_t       f_met;
   Float_t       f_tcmet;
   Float_t       f_pfmet;
   Float_t       f_mt;
   Float_t       f_tcmt;
   Float_t       f_pfmt;
   Float_t       f_pt;
   Float_t       f_eta;
   Float_t       f_phi;
   Int_t         f_charge;
   Float_t       f_weight;
   Int_t         f_nPFJets;
   Int_t         f_nJets;
   Float_t       f_event;
   Float_t       f_isIdWP70;
   Float_t       f_isIdWP80;
   Float_t       f_isIdWP85;
   Float_t       f_isIdWP90;
   Float_t       f_isIdWP95;
   Float_t       f_isIsoWP70;
   Float_t       f_isIsoWP80;
   Float_t       f_isIsoWP85;
   Float_t       f_isIsoWP90;
   Float_t       f_isIsoWP95;
   Float_t       f_isConvRejWP70;
   Float_t       f_isConvRejWP80;
   Float_t       f_isConvRejWP85;
   Float_t       f_isConvRejWP90;
   Float_t       f_isConvRejWP95;
   Float_t       f_isWP70;
   Float_t       f_isWP80;
   Float_t       f_isWP85;
   Float_t       f_isWP90;
   Float_t       f_isWP95;

   // Declaration of leaf types for Z T & P
   Double_t      ztap_N_sig_sw;
   Double_t      ztap_L_N_sig;
   Double_t      ztap_N_bkg_sw;
   Double_t      ztap_L_N_bkg;
   Double_t      ztap_eopout;
   Double_t      ztap_eop;
   Double_t      ztap_esc;
   Double_t      ztap_pin;
   Double_t      ztap_hoe;
   Double_t      ztap_deta;
   Double_t      ztap_dphi;
   Double_t      ztap_s9s25;
   Double_t      ztap_s1s9;
   Double_t      ztap_see;
   Double_t      ztap_fbrem;
   Double_t      ztap_missingHits;
   Double_t      ztap_convDcot;
   Double_t      ztap_convDist;
   Double_t      ztap_zmass;
   Double_t      ztap_charge;
   Double_t      ztap_eta;
   Double_t      ztap_pt;
   Double_t      ztap_iecal;
   Double_t      ztap_iptbin;
   Double_t      ztap_iclass;
   Double_t      ztap_weight;

   // List of branches
   TBranch        *b_nPU;
   TBranch        *b_N_sig_sw;   //!
   TBranch        *b_L_N_sig;   //!
   TBranch        *b_N_qcd_sw;   //!
   TBranch        *b_L_N_qcd;   //!
   TBranch        *b_trackerIso;   //!
   TBranch        *b_ecalJIso;   //!
   TBranch        *b_ecalGTIso;   //!
   TBranch        *b_hcalIso;   //!
   TBranch        *b_combinedIso;   //!
   TBranch        *b_classification;   //!
   TBranch        *b_deta;   //!
   TBranch        *b_dphi;   //!
   TBranch        *b_hoe;   //!
   TBranch        *b_see;   //!
   TBranch        *b_spp;   //!
   TBranch        *b_eop;   //!
   TBranch        *b_esc;   //!
   TBranch        *b_pin;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_fbrem;   //!
   TBranch        *b_nbrem;   //!
   TBranch        *b_met;   //!
   TBranch        *b_tcmet;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_mt;   //!
   TBranch        *b_tcmt;   //!
   TBranch        *b_pfmt;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_nPFJets;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_event;   //!
   TBranch        *b_isIdWP70;   //!
   TBranch        *b_isIdWP80;   //!
   TBranch        *b_isIdWP85;   //!
   TBranch        *b_isIdWP90;   //!
   TBranch        *b_isIdWP95;   //!
   TBranch        *b_isIsoWP70;   //!
   TBranch        *b_isIsoWP80;   //!
   TBranch        *b_isIsoWP85;   //!
   TBranch        *b_isIsoWP90;   //!
   TBranch        *b_isIsoWP95;   //!
   TBranch        *b_isConvRejWP70;   //!
   TBranch        *b_isConvRejWP80;   //!
   TBranch        *b_isConvRejWP85;   //!
   TBranch        *b_isConvRejWP90;   //!
   TBranch        *b_isConvRejWP95;   //!
   TBranch        *b_isWP70;   //!
   TBranch        *b_isWP80;   //!
   TBranch        *b_isWP85;   //!
   TBranch        *b_isWP90;   //!
   TBranch        *b_isWP95;   //!

   // List of branches for Z T & P
   TBranch        *b_N_bkg_sw;   //!
   TBranch        *b_L_N_bkg;   //!
   TBranch        *b_EoPout;   //!
   TBranch        *b_EoP;   //!
   TBranch        *b_Esc;   //!
   TBranch        *b_Pin;   //!
   TBranch        *b_HoE;   //!
   TBranch        *b_deltaEtaCorr;   //!
   TBranch        *b_deltaPhiCorr;   //!
   TBranch        *b_s9s25;   //!
   TBranch        *b_s1s9;   //!
   TBranch        *b_sigmaIEtaIEta;   //!
   TBranch        *b_fBrem;   //!
   TBranch        *b_missingHits;   //!
   TBranch        *b_convDcot;   //!
   TBranch        *b_convDist;   //!
   TBranch        *b_zmass;   //!
   TBranch        *b_iecal;   //!
   TBranch        *b_iptbin;   //!
   TBranch        *b_iclass;   //!

   sPlotsPdfsComparison();
   virtual ~sPlotsPdfsComparison();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *treel, int isMC=1, int isZTaP=0);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     bookHistosVariableBinning();
   virtual void     bookHistosFixedBinning();
   virtual void     bookFullHistos();
   virtual void     doSignalsPlots(bool what) { m_doSignal = what; }
   virtual float    likelihoodRatio(int isMc, ElectronLikelihood &lh);

 protected:

   bool m_isMC;
   bool m_doSignal;
   bool m_isDataZTaP;

   TH1F *etaEle;   
   TH1F *dPhiEle[3];
   TH1F *dEtaEle[3];
   TH1F *EoPEle[3];
   TH1F *OneOverEMinusOneOverPEle[3];
   TH1F *HoEEle[3];
   TH1F *sigmaIEtaIEtaEle[3];
   TH1F *fbremEle[3];
   TH1F *phiEle[3];
   TH1F *chargeEle[3];
   TH1F *lhEle[3];

  // ---------- monitoring histograms ------------

  /// Electrons: not splitted
  /// histo[ecalsubdet][ptbin]
  TH1F *dPhiUnsplitEle[3][2];
  TH1F *dEtaUnsplitEle[3][2];
  TH1F *EoPUnsplitEle[3][2];
  TH1F *OneOverEMinusOneOverPUnsplitEle[3][2];
  TH1F *HoEUnsplitEle[3][2];  
  TH1F *sigmaIEtaIEtaUnsplitEle[3][2];
  TH1F *sigmaIPhiIPhiUnsplitEle[3][2];
  TH1F *fBremUnsplitEle[3][2];
  TH1F *lhUnsplitEle[3][2];
  
  /// Electrons class-splitted
  /// histo[ecalsubdet][ptbin][class]
  TH1F *dPhiClassEle[3][2][2];
  TH1F *dEtaClassEle[3][2][2];
  TH1F *EoPClassEle[3][2][2];
  TH1F *OneOverEMinusOneOverPClassEle[3][2][2];
  TH1F *HoEClassEle[3][2][2];
  TH1F *sigmaIEtaIEtaClassEle[3][2][2];
  TH1F *sigmaIPhiIPhiClassEle[3][2][2];
  TH1F *fBremClassEle[3][2][2];
  TH1F *lhClassEle[3][2][2];

  // the likelihood algorithm
  ElectronLikelihood *LH;

};

#endif

#ifdef sPlotsPdfsComparison_cxx
sPlotsPdfsComparison::sPlotsPdfsComparison() {}

sPlotsPdfsComparison::~sPlotsPdfsComparison()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sPlotsPdfsComparison::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sPlotsPdfsComparison::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void sPlotsPdfsComparison::Init(TTree *tree, int isMC, int data_ZTaP)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   if(isMC) {
     std::cout << "Setting branches for MC tree" << std::endl;
     m_isMC = 1;
     fChain->SetBranchAddress("nPU", &f_nPU, &b_nPU);
     fChain->SetBranchAddress("trackerIso", &f_trackerIso, &b_trackerIso);
     fChain->SetBranchAddress("ecalJIso", &f_ecalJIso, &b_ecalJIso);
     fChain->SetBranchAddress("ecalGTIso", &f_ecalGTIso, &b_ecalGTIso);
     fChain->SetBranchAddress("hcalIso", &f_hcalIso, &b_hcalIso);
     fChain->SetBranchAddress("combinedIso", &f_combinedIso, &b_combinedIso);
     fChain->SetBranchAddress("classification", &f_classification, &b_classification);
     fChain->SetBranchAddress("deta", &f_deta, &b_deta);
     fChain->SetBranchAddress("dphi", &f_dphi, &b_dphi);
     fChain->SetBranchAddress("hoe", &f_hoe, &b_hoe);
     fChain->SetBranchAddress("see", &f_see, &b_see);
     fChain->SetBranchAddress("spp", &f_spp, &b_spp);
     fChain->SetBranchAddress("eop", &f_eop, &b_eop);
     fChain->SetBranchAddress("esc", &f_esc, &b_esc);
     fChain->SetBranchAddress("pin", &f_pin, &b_pin);
     fChain->SetBranchAddress("fbrem", &f_fbrem, &b_fbrem);
     fChain->SetBranchAddress("nbrem", &f_nbrem, &b_nbrem);
     fChain->SetBranchAddress("met", &f_met, &b_met);
     fChain->SetBranchAddress("tcmet", &f_tcmet, &b_tcmet);
     fChain->SetBranchAddress("pfmet", &f_pfmet, &b_pfmet);
     fChain->SetBranchAddress("mt", &f_mt, &b_mt);
     fChain->SetBranchAddress("tcmt", &f_tcmt, &b_tcmt);
     fChain->SetBranchAddress("pfmt", &f_pfmt, &b_pfmt);
     fChain->SetBranchAddress("pt1", &f_pt, &b_pt);
     fChain->SetBranchAddress("eta1", &f_eta, &b_eta);
     fChain->SetBranchAddress("phi1", &f_phi, &b_phi);
     fChain->SetBranchAddress("charge", &f_charge, &b_charge);
     fChain->SetBranchAddress("weight", &f_weight, &b_weight);
     fChain->SetBranchAddress("event", &f_event, &b_event);
     fChain->SetBranchAddress("nPFJetsHi", &f_nPFJets, &b_nPFJets);
     fChain->SetBranchAddress("nJetsHi", &f_nJets, &b_nJets);
     fChain->SetBranchAddress("isIdWP70", &f_isIdWP70, &b_isIdWP70);
     fChain->SetBranchAddress("isIdWP80", &f_isIdWP80, &b_isIdWP80);
     fChain->SetBranchAddress("isIdWP85", &f_isIdWP85, &b_isIdWP85);
     fChain->SetBranchAddress("isIdWP90", &f_isIdWP90, &b_isIdWP90);
     fChain->SetBranchAddress("isIdWP95", &f_isIdWP95, &b_isIdWP95);
     fChain->SetBranchAddress("isIsoWP70", &f_isIsoWP70, &b_isIsoWP70);
     fChain->SetBranchAddress("isIsoWP80", &f_isIsoWP80, &b_isIsoWP80);
     fChain->SetBranchAddress("isIsoWP85", &f_isIsoWP85, &b_isIsoWP85);
     fChain->SetBranchAddress("isIsoWP90", &f_isIsoWP90, &b_isIsoWP90);
     fChain->SetBranchAddress("isIsoWP95", &f_isIsoWP95, &b_isIsoWP95);
     fChain->SetBranchAddress("isConvRejWP70", &f_isConvRejWP70, &b_isConvRejWP70);
     fChain->SetBranchAddress("isConvRejWP80", &f_isConvRejWP80, &b_isConvRejWP80);
     fChain->SetBranchAddress("isConvRejWP85", &f_isConvRejWP85, &b_isConvRejWP85);
     fChain->SetBranchAddress("isConvRejWP90", &f_isConvRejWP90, &b_isConvRejWP90);
     fChain->SetBranchAddress("isConvRejWP95", &f_isConvRejWP95, &b_isConvRejWP95);
     fChain->SetBranchAddress("isWP70", &f_isWP70, &b_isWP70);
     fChain->SetBranchAddress("isWP80", &f_isWP80, &b_isWP80);
     fChain->SetBranchAddress("isWP85", &f_isWP85, &b_isWP85);
     fChain->SetBranchAddress("isWP90", &f_isWP90, &b_isWP90);
     fChain->SetBranchAddress("isWP95", &f_isWP95, &b_isWP95);
   } else {
     m_isMC = 0;
     fChain->SetBranchAddress("N_sig_sw", &N_sig_sw, &b_N_sig_sw);
     fChain->SetBranchAddress("L_N_sig", &L_N_sig, &b_L_N_sig);
     fChain->SetBranchAddress("N_qcd_sw", &N_qcd_sw, &b_N_qcd_sw);
     fChain->SetBranchAddress("L_N_qcd", &L_N_qcd, &b_L_N_qcd);
     fChain->SetBranchAddress("trackerIso", &trackerIso, &b_trackerIso);
     fChain->SetBranchAddress("ecalJIso", &ecalJIso, &b_ecalJIso);
     fChain->SetBranchAddress("ecalGTIso", &ecalGTIso, &b_ecalGTIso);
     fChain->SetBranchAddress("hcalIso", &hcalIso, &b_hcalIso);
     fChain->SetBranchAddress("combinedIso", &combinedIso, &b_combinedIso);
     fChain->SetBranchAddress("classification", &classification, &b_classification);
     fChain->SetBranchAddress("deta", &deta, &b_deta);
     fChain->SetBranchAddress("dphi", &dphi, &b_dphi);
     fChain->SetBranchAddress("hoe", &hoe, &b_hoe);
     fChain->SetBranchAddress("see", &see, &b_see);
     fChain->SetBranchAddress("spp", &spp, &b_spp);
     fChain->SetBranchAddress("eop", &eop, &b_eop);
     fChain->SetBranchAddress("esc", &esc, &b_esc);
     fChain->SetBranchAddress("pin", &pin, &b_pin);
     fChain->SetBranchAddress("fbrem", &fbrem, &b_fbrem);
     fChain->SetBranchAddress("nbrem", &nbrem, &b_nbrem);
     fChain->SetBranchAddress("met", &met, &b_met);
     fChain->SetBranchAddress("tcmet", &tcmet, &b_tcmet);
     fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
     fChain->SetBranchAddress("mt", &mt, &b_mt);
     fChain->SetBranchAddress("tcmt", &tcmt, &b_tcmt);
     fChain->SetBranchAddress("pfmt", &pfmt, &b_pfmt);
     fChain->SetBranchAddress("pt", &pt, &b_pt);
     fChain->SetBranchAddress("eta", &eta, &b_eta);
     fChain->SetBranchAddress("phi", &phi, &b_phi);
     fChain->SetBranchAddress("charge", &charge, &b_charge);
     fChain->SetBranchAddress("weight", &weight, &b_weight);
     fChain->SetBranchAddress("event", &event, &b_event);
     fChain->SetBranchAddress("nPFJetsHi", &nPFJets, &b_nPFJets);
     fChain->SetBranchAddress("nJetsHi", &nJets, &b_nJets);
     fChain->SetBranchAddress("isIdWP70", &isIdWP70, &b_isIdWP70);
     fChain->SetBranchAddress("isIdWP80", &isIdWP80, &b_isIdWP80);
     fChain->SetBranchAddress("isIdWP85", &isIdWP85, &b_isIdWP85);
     fChain->SetBranchAddress("isIdWP90", &isIdWP90, &b_isIdWP90);
     fChain->SetBranchAddress("isIdWP95", &isIdWP95, &b_isIdWP95);
     fChain->SetBranchAddress("isIsoWP70", &isIsoWP70, &b_isIsoWP70);
     fChain->SetBranchAddress("isIsoWP80", &isIsoWP80, &b_isIsoWP80);
     fChain->SetBranchAddress("isIsoWP85", &isIsoWP85, &b_isIsoWP85);
     fChain->SetBranchAddress("isIsoWP90", &isIsoWP90, &b_isIsoWP90);
     fChain->SetBranchAddress("isIsoWP95", &isIsoWP95, &b_isIsoWP95);
     fChain->SetBranchAddress("isConvRejWP70", &isConvRejWP70, &b_isConvRejWP70);
     fChain->SetBranchAddress("isConvRejWP80", &isConvRejWP80, &b_isConvRejWP80);
     fChain->SetBranchAddress("isConvRejWP85", &isConvRejWP85, &b_isConvRejWP85);
     fChain->SetBranchAddress("isConvRejWP90", &isConvRejWP90, &b_isConvRejWP90);
     fChain->SetBranchAddress("isConvRejWP95", &isConvRejWP95, &b_isConvRejWP95);
     fChain->SetBranchAddress("isWP70", &isWP70, &b_isWP70);
     fChain->SetBranchAddress("isWP80", &isWP80, &b_isWP80);
     fChain->SetBranchAddress("isWP85", &isWP85, &b_isWP85);
     fChain->SetBranchAddress("isWP90", &isWP90, &b_isWP90);
     fChain->SetBranchAddress("isWP95", &isWP95, &b_isWP95);
   }

   if(data_ZTaP) {
     m_isDataZTaP = 1;
     fChain->SetBranchAddress("N_sig_sw", &ztap_N_sig_sw, &b_N_sig_sw);
     fChain->SetBranchAddress("L_N_sig", &ztap_L_N_sig, &b_L_N_sig);
     fChain->SetBranchAddress("N_bkg_sw", &ztap_N_bkg_sw, &b_N_bkg_sw);
     fChain->SetBranchAddress("L_N_bkg", &ztap_L_N_bkg, &b_L_N_bkg);
     fChain->SetBranchAddress("EoPout", &ztap_eopout, &b_EoPout);
     fChain->SetBranchAddress("EoP", &ztap_eop, &b_EoP);
     fChain->SetBranchAddress("esc", &ztap_esc, &b_Esc);
     fChain->SetBranchAddress("pin", &ztap_pin, &b_Pin);
     fChain->SetBranchAddress("HoE", &ztap_hoe, &b_HoE);
     fChain->SetBranchAddress("deltaEtaCorr", &ztap_deta, &b_deltaEtaCorr);
     fChain->SetBranchAddress("deltaPhiCorr", &ztap_dphi, &b_deltaPhiCorr);
     fChain->SetBranchAddress("s9s25", &ztap_s9s25, &b_s9s25);
     fChain->SetBranchAddress("s1s9", &ztap_s1s9, &b_s1s9);
     fChain->SetBranchAddress("sigmaIEtaIEta", &ztap_see, &b_sigmaIEtaIEta);
     fChain->SetBranchAddress("fBrem", &ztap_fbrem, &b_fBrem);
     fChain->SetBranchAddress("missingHits", &ztap_missingHits, &b_missingHits);
     fChain->SetBranchAddress("convDcot", &ztap_convDcot, &b_convDcot);
     fChain->SetBranchAddress("convDist", &ztap_convDist, &b_convDist);
     fChain->SetBranchAddress("zmass", &ztap_zmass, &b_zmass);
     fChain->SetBranchAddress("charge", &ztap_charge, &b_charge);
     fChain->SetBranchAddress("eta", &ztap_eta, &b_eta);
     fChain->SetBranchAddress("pt", &ztap_pt, &b_pt);
     fChain->SetBranchAddress("iecal", &ztap_iecal, &b_iecal);
     fChain->SetBranchAddress("iptbin", &ztap_iptbin, &b_iptbin);
     fChain->SetBranchAddress("iclass", &ztap_iclass, &b_iclass);
     fChain->SetBranchAddress("weight", &ztap_weight, &b_weight);
   } else {m_isDataZTaP=0;}
   Notify();
}

Bool_t sPlotsPdfsComparison::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sPlotsPdfsComparison::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sPlotsPdfsComparison::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sPlotsPdfsComparison_cxx
