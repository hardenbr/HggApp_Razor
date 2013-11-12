//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct  9 11:32:05 2009 by ROOT version 5.22/00a
// from TTree dataset/dataset with sWeights
// found on file: /cmsrm/pc21/emanuele/data/Likelihood3.1.4/sPlotsTrees/sPlotsZee_tree.root
//////////////////////////////////////////////////////////

#ifndef sPlotsPdfsComparison_h
#define sPlotsPdfsComparison_h

#include <TROOT.h>
#include <TChain.h>
#include <TH1F.h>
#include <TFile.h>

class sPlotsPdfsComparison {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types: NORMAL TREE! contains floats, ints, etc...
//    Double_t        N_sig_sw;
//    Double_t        L_N_sig;
//    Double_t        N_bkg_sw;
//    Double_t        L_N_bkg;
//    Float_t        EoPout;
//    Float_t        EoP;
//    Float_t        HoE;
//    Float_t        deltaEta;
//    Float_t        deltaPhi;
//    Float_t        s9s25;
//    Float_t        s1s9;
//    Float_t        sigmaIEtaIEta;
//    Float_t        zmass;
//    Float_t        charge;
//    Float_t        eta;
//    Float_t        pt;
//    Int_t        iecal;
//    Int_t        iptbin;
//    Int_t        iclass;
//    Float_t        weight;

   // Declaration of leaf types: TREE FROM DATASET! contains only doubles (Roofit grrrrrrrrr....)
   Double_t        N_sig_sw;
   Double_t        L_N_sig;
   Double_t        N_bkg_sw;
   Double_t        L_N_bkg;
   Double_t        EoPout;
   Double_t        EoP;
   Double_t        HoE;
   Double_t        deltaEta;
   Double_t        deltaPhi;
   Double_t        s9s25;
   Double_t        s1s9;
   Double_t        sigmaIEtaIEta;
   Double_t        zmass;
   Double_t        charge;
   Double_t        eta;
   Double_t        pt;
   Double_t        iecal;
   Double_t        iptbin;
   Double_t        iclass;
   Double_t        weight;

   // List of branches
   TBranch        *b_N_sig_sw;   //!
   TBranch        *b_L_N_sig;   //!
   TBranch        *b_N_bkg_sw;   //!
   TBranch        *b_L_N_bkg;   //!
   TBranch        *b_EoPout;   //!
   TBranch        *b_EoP;   //!
   TBranch        *b_HoE;   //!
   TBranch        *b_deltaEta;   //!
   TBranch        *b_deltaPhi;   //!
   TBranch        *b_s9s25;   //!
   TBranch        *b_s1s9;   //!
   TBranch        *b_sigmaIEtaIEta;   //!
   TBranch        *b_zmass;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_iecal;   //!
   TBranch        *b_iptbin;   //!
   TBranch        *b_iclass;   //!
   TBranch        *b_weight;   //!

   sPlotsPdfsComparison(TTree *tree=0);
   virtual ~sPlotsPdfsComparison();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     RunOverMC(bool what) { m_isMC = what; }
   virtual void     bookHistos();

 protected:
   bool m_isMC;

  /// Electrons class-splitted
  /// histo[ecalsubdet][ptbin][class]
  TH1F *dPhiClassEle[2][2][2];
  TH1F *dEtaClassEle[2][2][2];
  TH1F *EoPoutClassEle[2][2][2];
  TH1F *EoPClassEle[2][2][2];
  TH1F *HoEClassEle[2][2][2];
  TH1F *sigmaIEtaIEtaClassEle[2][2][2];
  TH1F *s1s9ClassEle[2][2][2];
  TH1F *s9s25ClassEle[2][2][2];

};

#endif

#ifdef sPlotsPdfsComparison_cxx
sPlotsPdfsComparison::sPlotsPdfsComparison(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cmsrm/pc21/emanuele/data/Likelihood3.1.4/sPlotsTrees/sPlotsZee_tree.root");
      if (!f) {
         f = new TFile("/cmsrm/pc21/emanuele/data/Likelihood3.1.4/sPlotsTrees/sPlotsZee_tree.root");
      }
      tree = (TTree*)gDirectory->Get("T1");

   }
   Init(tree);
}

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

void sPlotsPdfsComparison::Init(TTree *tree)
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

   fChain->SetBranchAddress("N_sig_sw", &N_sig_sw, &b_N_sig_sw);
   fChain->SetBranchAddress("L_N_sig", &L_N_sig, &b_L_N_sig);
   fChain->SetBranchAddress("N_bkg_sw", &N_bkg_sw, &b_N_bkg_sw);
   fChain->SetBranchAddress("L_N_bkg", &L_N_bkg, &b_L_N_bkg);
   fChain->SetBranchAddress("EoPout", &EoPout, &b_EoPout);
   fChain->SetBranchAddress("EoP", &EoP, &b_EoP);
   fChain->SetBranchAddress("HoE", &HoE, &b_HoE);
   fChain->SetBranchAddress("deltaEta", &deltaEta, &b_deltaEta);
   fChain->SetBranchAddress("deltaPhi", &deltaPhi, &b_deltaPhi);
   fChain->SetBranchAddress("s9s25", &s9s25, &b_s9s25);
   fChain->SetBranchAddress("s1s9", &s1s9, &b_s1s9);
   fChain->SetBranchAddress("sigmaIEtaIEta", &sigmaIEtaIEta, &b_sigmaIEtaIEta);
   fChain->SetBranchAddress("zmass", &zmass, &b_zmass);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("eta", &eta, &b_eta);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("iecal", &iecal, &b_iecal);
   fChain->SetBranchAddress("iptbin", &iptbin, &b_iptbin);
   fChain->SetBranchAddress("iclass", &iclass, &b_iclass);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
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
