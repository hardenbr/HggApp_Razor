// --- usage: ./CompareClasses EleEfficiencyEta.root EleMisidEta.root EleEfficiencyPt.root EleMisidPt.root
// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>

// Offline analysis includes
#include "CommonTools/include/EfficiencyEvaluator.hh"

using namespace std;

void drawEta(const char* filenameSig, const char* filenameBkg);
void drawPt(const char* filenameSig, const char* filenameBkg);

int main(int argc, char* argv[]) {

  char inputFileNameEtaSig[150], inputFileNameEtaBkg[150];
  char inputFileNamePtSig[150], inputFileNamePtBkg[150];
  if ( argc < 5 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: compareMisid inputFileEtaSig.root inputFileEtaBkg.root inputFilePtSig.root inputFilePtBkg.root" << std::endl;
    return 1;
  }
  strcpy(inputFileNameEtaSig,argv[1]);
  strcpy(inputFileNameEtaBkg,argv[2]);
  strcpy(inputFileNamePtSig,argv[3]);
  strcpy(inputFileNamePtBkg,argv[4]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  drawEta(inputFileNameEtaSig, inputFileNameEtaBkg);
  drawPt(inputFileNamePtSig, inputFileNamePtBkg);

}

void drawEta(const char* filenameSig, const char* filenameBkg) {

  TFile *efficiencyFileEtaSig = TFile::Open(filenameSig);

  TH1F *RecoEtaSig = (TH1F*)efficiencyFileEtaSig->Get("RecoEta");
  TH1F *GoldenEtaSig = (TH1F*)efficiencyFileEtaSig->Get("GoldenEta");
  TH1F *BigBremEtaSig = (TH1F*)efficiencyFileEtaSig->Get("BigBremEta");
  TH1F *NarrowEtaSig = (TH1F*)efficiencyFileEtaSig->Get("NarrowEta");
  TH1F *ShoweringEtaSig = (TH1F*)efficiencyFileEtaSig->Get("ShoweringEta");

  TH1F *NonShoweringEtaSig = (TH1F*)GoldenEtaSig->Clone("NonShoweringEtaSig");
  NonShoweringEtaSig->Add(GoldenEtaSig);
  NonShoweringEtaSig->Add(BigBremEtaSig);
  NonShoweringEtaSig->Add(NarrowEtaSig);
  
  EfficiencyEvaluator ElectronClassesEtaSig("fakesclassessig-eta.root");
  ElectronClassesEtaSig.AddNumerator(RecoEtaSig);
  ElectronClassesEtaSig.AddNumerator(ShoweringEtaSig);
  ElectronClassesEtaSig.SetDenominator(RecoEtaSig);
  ElectronClassesEtaSig.ComputeEfficiencies();

  vector<TH1F*> fractionSig = ElectronClassesEtaSig.GetCumulativeEfficiencies();

  TFile *efficiencyFileEtaBkg = TFile::Open(filenameBkg);

  TH1F *RecoEtaBkg = (TH1F*)efficiencyFileEtaBkg->Get("RecoEta");
  TH1F *GoldenEtaBkg = (TH1F*)efficiencyFileEtaBkg->Get("GoldenEta");
  TH1F *BigBremEtaBkg = (TH1F*)efficiencyFileEtaBkg->Get("BigBremEta");
  TH1F *NarrowEtaBkg = (TH1F*)efficiencyFileEtaBkg->Get("NarrowEta");
  TH1F *ShoweringEtaBkg = (TH1F*)efficiencyFileEtaBkg->Get("ShoweringEta");

  TH1F *NonShoweringEtaBkg = (TH1F*)GoldenEtaBkg->Clone("NonShoweringEtaBkg");
  NonShoweringEtaBkg->Add(GoldenEtaBkg);
  NonShoweringEtaBkg->Add(BigBremEtaBkg);
  NonShoweringEtaBkg->Add(NarrowEtaBkg);
  
  EfficiencyEvaluator ElectronClassesEtaBkg("fakesclassessig-eta.root");
  ElectronClassesEtaBkg.AddNumerator(RecoEtaBkg);
  ElectronClassesEtaBkg.AddNumerator(ShoweringEtaBkg);
  ElectronClassesEtaBkg.SetDenominator(RecoEtaBkg);
  ElectronClassesEtaBkg.ComputeEfficiencies();

  vector<TH1F*> fractionBkg = ElectronClassesEtaBkg.GetCumulativeEfficiencies();

  TCanvas c1;
  c1.SetLogy(0);
  fractionSig[1]->SetLineColor(2);
  fractionBkg[1]->SetLineColor(4);
  
  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(fractionSig[1],"electrons");
  leg->AddEntry(fractionBkg[1],"jets");

  fractionSig[1]->SetMinimum(0.0);
  fractionSig[1]->SetMaximum(1.4);
  fractionSig[1]->SetTitle("");
  fractionSig[1]->GetXaxis()->SetTitle("electron #eta");
  fractionSig[1]->GetYaxis()->SetTitle("showering+cracks fraction");
  fractionSig[1]->Draw("hist");
  fractionBkg[1]->Draw("same hist");

  leg->Draw();

  c1.SaveAs("classfraction-eta.eps");
  c1.SaveAs("classfraction-eta.root");

}


void drawPt(const char* filenameSig, const char* filenameBkg) {

  TFile *efficiencyFilePtSig = TFile::Open(filenameSig);

  TH1F *RecoPtSig = (TH1F*)efficiencyFilePtSig->Get("RecoPt");
  TH1F *GoldenPtSig = (TH1F*)efficiencyFilePtSig->Get("GoldenPt");
  TH1F *BigBremPtSig = (TH1F*)efficiencyFilePtSig->Get("BigBremPt");
  TH1F *NarrowPtSig = (TH1F*)efficiencyFilePtSig->Get("NarrowPt");
  TH1F *ShoweringPtSig = (TH1F*)efficiencyFilePtSig->Get("ShoweringPt");

  TH1F *NonShoweringPtSig = (TH1F*)GoldenPtSig->Clone("NonShoweringPtSig");
  NonShoweringPtSig->Add(GoldenPtSig);
  NonShoweringPtSig->Add(BigBremPtSig);
  NonShoweringPtSig->Add(NarrowPtSig);
  
  EfficiencyEvaluator ElectronClassesPtSig("fakesclassessig-pt.root");
  ElectronClassesPtSig.AddNumerator(RecoPtSig);
  ElectronClassesPtSig.AddNumerator(ShoweringPtSig);
  ElectronClassesPtSig.SetDenominator(RecoPtSig);
  ElectronClassesPtSig.ComputeEfficiencies();

  vector<TH1F*> fractionSig = ElectronClassesPtSig.GetCumulativeEfficiencies();

  TFile *efficiencyFilePtBkg = TFile::Open(filenameBkg);

  TH1F *RecoPtBkg = (TH1F*)efficiencyFilePtBkg->Get("RecoPt");
  TH1F *GoldenPtBkg = (TH1F*)efficiencyFilePtBkg->Get("GoldenPt");
  TH1F *BigBremPtBkg = (TH1F*)efficiencyFilePtBkg->Get("BigBremPt");
  TH1F *NarrowPtBkg = (TH1F*)efficiencyFilePtBkg->Get("NarrowPt");
  TH1F *ShoweringPtBkg = (TH1F*)efficiencyFilePtBkg->Get("ShoweringPt");

  TH1F *NonShoweringPtBkg = (TH1F*)GoldenPtBkg->Clone("NonShoweringPtBkg");
  NonShoweringPtBkg->Add(GoldenPtBkg);
  NonShoweringPtBkg->Add(BigBremPtBkg);
  NonShoweringPtBkg->Add(NarrowPtBkg);
  
  EfficiencyEvaluator ElectronClassesPtBkg("fakesclassessig-eta.root");
  ElectronClassesPtBkg.AddNumerator(RecoPtBkg);
  ElectronClassesPtBkg.AddNumerator(ShoweringPtBkg);
  ElectronClassesPtBkg.SetDenominator(RecoPtBkg);
  ElectronClassesPtBkg.ComputeEfficiencies();

  vector<TH1F*> fractionBkg = ElectronClassesPtBkg.GetCumulativeEfficiencies();

  TCanvas c1;
  c1.SetLogy(0);
  fractionSig[1]->SetLineColor(2);
  fractionBkg[1]->SetLineColor(4);
  
  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(fractionSig[1],"electrons");
  leg->AddEntry(fractionBkg[1],"jets");

  fractionSig[1]->SetMinimum(0.0);
  fractionSig[1]->SetMaximum(1.4);
  fractionSig[1]->SetTitle("");
  fractionSig[1]->GetXaxis()->SetTitle("electron p_{T} [GeV]");
  fractionSig[1]->GetYaxis()->SetTitle("showering+cracks fraction");
  fractionSig[1]->Draw("hist");
  fractionBkg[1]->Draw("same hist");

  leg->Draw();

  c1.SaveAs("classfraction-pt.eps");
  c1.SaveAs("classfraction-pt.root");

}
