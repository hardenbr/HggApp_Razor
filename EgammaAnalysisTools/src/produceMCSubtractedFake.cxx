// ROOT includes
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include "TStyle.h"
#include "TROOT.h"

// C++ includes
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Offline analysis includes                
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "TString.h"
#include <string>

using namespace std;

void computeAverage(vector<TH1F*> efficiencyLH);

// Data histos in input must not be weighted
//   MC histos in input must have a weight corresponding to 1/pb 

int main(int argc, char* argv[]) {
  
  if ( argc < 2 ){
    std::cout << "missing argument! insert lumi" << std::endl;
    return 1;
  }
  float lumi = atof(argv[1]);
  cout << "rescale for lumi = " << lumi << endl;

  // input files
  TFile *filePT_EB[3];
  filePT_EB[0] = new TFile("results_data/EleMisidPtBarrel2.root");
  // filePT_EB[0] = new TFile("macro/fakeTagAndProbe/newWeightedHistos_1p_mergedSel_QCD/mergedResultPtBarrel2QCD.root");
  filePT_EB[1] = new TFile("macro/fakeTagAndProbe/newWeightedHistos_1p_mergedSel_or2HLT_ET35_2011AandB/mergedResultPtBarrel2Wjets.root");
  filePT_EB[2] = new TFile("macro/fakeTagAndProbe/newWeightedHistos_1p_mergedSel_or2HLT_ET35_2011AandB/mergedResultPtBarrel2DY.root");

  TFile *filePT_EE[3];
  filePT_EE[0] = new TFile("results_data/EleMisidPtEndcap2.root");
  // filePT_EE[0] = new TFile("macro/fakeTagAndProbe/newWeightedHistos_1p_mergedSel_QCD/mergedResultPtEndcap2QCD.root");
  filePT_EE[1] = new TFile("macro/fakeTagAndProbe/newWeightedHistos_1p_mergedSel_or2HLT_ET35_2011AandB/mergedResultPtEndcap2Wjets.root");
  filePT_EE[2] = new TFile("macro/fakeTagAndProbe/newWeightedHistos_1p_mergedSel_or2HLT_ET35_2011AandB/mergedResultPtEndcap2DY.root");

  // taking histos before EWK subtraction 
  TH1F *LHIdLoosePt_PRE_EB[3], *FakeableJetsPt_PRE_EB[3];
  TH1F *LHIdLoosePt_PRE_EE[3], *FakeableJetsPt_PRE_EE[3];
  for (int ii=0; ii<3; ii++) { 
    cout << ii << endl;
    LHIdLoosePt_PRE_EB[ii] = (TH1F*)filePT_EB[ii]->Get("BdtIdWP80SmurfPtBarrel2");
    FakeableJetsPt_PRE_EB[ii] = (TH1F*)filePT_EB[ii]->Get("RecoPtBarrel2");
    // 
    cout << "presi i barrel" << endl;
    LHIdLoosePt_PRE_EE[ii] = (TH1F*)filePT_EE[ii]->Get("BdtIdWP80SmurfPtEndcap2");    
    cout << "done" << endl;
    FakeableJetsPt_PRE_EE[ii] = (TH1F*)filePT_EE[ii]->Get("RecoPtEndcap2");
  }
  cout << "histos" << endl;

  // subtracting EWK from denominator
  TH1F *FakeableJetsPt_EB = (TH1F*)FakeableJetsPt_PRE_EB[0] -> Clone("FakeableJetsPt_EB");
  for (int ibin=0; ibin<10; ibin++) {
    float dataBinContent = FakeableJetsPt_PRE_EB[0]->GetBinContent(ibin);
    float WBinContent    = lumi * FakeableJetsPt_PRE_EB[1]->GetBinContent(ibin);
    float ZBinContent    = lumi * FakeableJetsPt_PRE_EB[2]->GetBinContent(ibin);
    float subtracted     = dataBinContent - WBinContent - ZBinContent;
    FakeableJetsPt_EB    -> SetBinContent( ibin, subtracted);
  }

  TH1F *FakeableJetsPt_EE = (TH1F*)FakeableJetsPt_PRE_EE[0] -> Clone("FakeableJetsPt_EE");
  for (int ibin=0; ibin<10; ibin++) {
    float dataBinContent = FakeableJetsPt_PRE_EE[0]->GetBinContent(ibin);
    float WBinContent    = lumi * FakeableJetsPt_PRE_EE[1]->GetBinContent(ibin);
    float ZBinContent    = lumi * FakeableJetsPt_PRE_EE[2]->GetBinContent(ibin);
    float subtracted     = dataBinContent - WBinContent - ZBinContent;
    FakeableJetsPt_EE    -> SetBinContent( ibin, subtracted);
  }

  // subtracting EWK from numerator
  TH1F *LHIdLoosePt_EB = (TH1F*)LHIdLoosePt_PRE_EB[0] -> Clone("BdtIdWP80SmurfPt_EB");  
  for (int ibin=0; ibin<10; ibin++) {
    float dataBinContent = LHIdLoosePt_PRE_EB[0]->GetBinContent(ibin);
    float WBinContent    = lumi * LHIdLoosePt_PRE_EB[1]->GetBinContent(ibin);
    float ZBinContent    = lumi * LHIdLoosePt_PRE_EB[2]->GetBinContent(ibin);
    float subtracted     = dataBinContent - WBinContent - ZBinContent;
    LHIdLoosePt_EB       -> SetBinContent( ibin, subtracted);
  }

  TH1F *LHIdLoosePt_EE = (TH1F*)LHIdLoosePt_PRE_EE[0] -> Clone("BdtIdWP80SmurfPt_EE");
  for (int ibin=0; ibin<10; ibin++) {
    float dataBinContent = LHIdLoosePt_PRE_EE[0]->GetBinContent(ibin);
    float WBinContent    = lumi * LHIdLoosePt_PRE_EE[1]->GetBinContent(ibin);
    float ZBinContent    = lumi * LHIdLoosePt_PRE_EE[2]->GetBinContent(ibin);
    float subtracted     = dataBinContent - WBinContent - ZBinContent;
    LHIdLoosePt_EE       -> SetBinContent( ibin, subtracted);
  }

  cout << "before any subtraction" << endl;
  for (int ibin=0; ibin<10; ibin++) cout << "data: bin " << ibin << ", content = " << LHIdLoosePt_PRE_EB[0]->GetBinContent(ibin) << endl;
  cout << endl;
  for (int ibin=0; ibin<10; ibin++) cout << "Wj:   bin " << ibin << ", content = " << LHIdLoosePt_PRE_EB[1]->GetBinContent(ibin) << endl;
  cout << endl;
  for (int ibin=0; ibin<10; ibin++) cout << "DY:   bin " << ibin << ", content = " << LHIdLoosePt_PRE_EB[2]->GetBinContent(ibin) << endl;
  cout << endl;
  cout << "after the subtraction" << endl;
  for (int ibin=0; ibin<10; ibin++) cout << "bin " << ibin << ", content = " << LHIdLoosePt_EB->GetBinContent(ibin) << endl;
  cout << endl;
  
  // now computing the fake rate - before and after subtraction
  EfficiencyEvaluator ElectronEfficiencyPtLH_EB("elefake-pt-barrel1.root");
  ElectronEfficiencyPtLH_EB.AddNumerator(LHIdLoosePt_EB);
  ElectronEfficiencyPtLH_EB.SetDenominator(FakeableJetsPt_EB);
  ElectronEfficiencyPtLH_EB.ComputeEfficiencies();
  vector<TH1F*> efficiencyLH_EB = ElectronEfficiencyPtLH_EB.GetCumulativeEfficiencies();
  
  EfficiencyEvaluator ElectronEfficiencyPtLH_EE("elefake-pt-endcap1.root");
  ElectronEfficiencyPtLH_EE.AddNumerator(LHIdLoosePt_EE);
  ElectronEfficiencyPtLH_EE.SetDenominator(FakeableJetsPt_EE);
  ElectronEfficiencyPtLH_EE.ComputeEfficiencies();
  vector<TH1F*> efficiencyLH_EE = ElectronEfficiencyPtLH_EE.GetCumulativeEfficiencies();

  cout << "barrel" << endl;
  computeAverage(efficiencyLH_EB);

  cout << endl;
  cout << "endcap" << endl;
  computeAverage(efficiencyLH_EE);

  // saving histos
  TFile file("myFR.root","RECREATE");
  efficiencyLH_EB[0]->Write();
  efficiencyLH_EE[0]->Write();
}


void computeAverage(vector<TH1F*> efficiencyLH) {

  int binLHA = efficiencyLH[0] -> GetNbinsX();
  int binLH  = binLHA+1;

  float contentLH[binLH], errorLH[binLH];

  for (int ibin=0; ibin<=binLH; ibin++) {
    contentLH[ibin] = efficiencyLH[0]->GetBinContent(ibin);
    if ( (efficiencyLH[0]->GetBinContent(ibin)) > 0 )
      errorLH[ibin] = efficiencyLH[0]->GetBinError(ibin);
    else
      errorLH[ibin] = 0.;
  }
  
  for (int ibin=0; ibin<=binLH; ibin++) cout << "bin = " << ibin << ", value = " << contentLH[ibin] << " +- " << errorLH[ibin] << endl;
}


