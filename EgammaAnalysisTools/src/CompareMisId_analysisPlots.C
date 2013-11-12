
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
#include "TString.h"
#include <string>

using namespace std;

std::vector<TString> EgammaLHBasedIDWPs;

void drawEta();
void drawPt();

int main(int argc, char* argv[]) {
  
  TString IDpart;
  if (argv[1]) {
    char idpar[100];
    strcpy(idpar,argv[1]);
    IDpart=TString(idpar);
  }
  else
    IDpart="";
  
  EgammaLHBasedIDWPs.push_back(TString("LHLoose")+TString(IDpart));  

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  
  drawEta();
  drawPt();
}


void drawEta() {
  
  TFile *efficiencyFileEta[6];
  // to compare HLT paths             
  efficiencyFileEta[0] = TFile::Open("results_data/allData_fakeTeP_may1/Data7TeV/DoubleElectron/4/EleMisidEta.root");
  efficiencyFileEta[1] = TFile::Open("results_data/allData_onlyEle8HLT_fakeTeP_may1/Data7TeV/DoubleElectron/4/EleMisidEta.root");
  efficiencyFileEta[2] = TFile::Open("results_data/allData_onlyEle17HLT_fakeTeP_may1/Data7TeV/DoubleElectron/4/EleMisidEta.root");
  // to compare data with mc
  efficiencyFileEta[3] = TFile::Open("macro/fakeTagAndProbe/lastMC_weightedHistos/mergedResultEtaQCD.root");
  efficiencyFileEta[4] = TFile::Open("macro/fakeTagAndProbe/lastMC_weightedHistos/mergedResultEtaWjets.root");
  efficiencyFileEta[5] = TFile::Open("macro/fakeTagAndProbe/lastMC_weightedHistos/mergedResultEtaDY.root");
  
  TH1F *RecoEta[6];
  for (int ii=0; ii<6; ii++) {
    RecoEta[ii] = (TH1F*)efficiencyFileEta[ii]->Get("FakeableJetsEta");
    RecoEta[ii]->Rebin(6);
  }

  std::vector<TH1F*> LHIdEta[6];
  for (int ii=0; ii<6; ii++) {
    for (int i=0;i<EgammaLHBasedIDWPs.size();++i) {
      TH1F* aHisto = (TH1F*)efficiencyFileEta[ii]->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta");
      aHisto->Rebin(6);
      LHIdEta[ii].push_back(aHisto);
    }
  }

  EfficiencyEvaluator *ElectronEfficiencyEta[6];
  ElectronEfficiencyEta[0] = new EfficiencyEvaluator("elefake-eta0.root");
  ElectronEfficiencyEta[1] = new EfficiencyEvaluator("elefake-eta1.root");
  ElectronEfficiencyEta[2] = new EfficiencyEvaluator("elefake-eta2.root");
  ElectronEfficiencyEta[3] = new EfficiencyEvaluator("elefake-eta3.root");
  ElectronEfficiencyEta[4] = new EfficiencyEvaluator("elefake-eta4.root");
  ElectronEfficiencyEta[5] = new EfficiencyEvaluator("elefake-eta5.root");

  for (int ii=0; ii<6; ii++) {
    for (int i=0;i<LHIdEta[ii].size();++i) ElectronEfficiencyEta[ii]->AddNumerator(LHIdEta[ii][i]);
    ElectronEfficiencyEta[ii]->SetDenominator(RecoEta[ii]);
    ElectronEfficiencyEta[ii]->ComputeEfficiencies();
  }

  vector<TH1F*> efficiency[6];
  for (int ii=0; ii<6; ii++) efficiency[ii] = ElectronEfficiencyEta[ii]->GetCumulativeEfficiencies();

  for (int tightness=0;tightness<LHIdEta[0].size();++tightness) {
    
    for (int ii=0; ii<6; ii++) {

      int theColor = ii+1;
      if (ii<3) theColor = ii+1;
      if (ii>2) theColor = ii-1;

      efficiency[ii][0+tightness]->SetLineColor(theColor);
      efficiency[ii][0+tightness]->SetMarkerColor(theColor);
      efficiency[ii][0+tightness]->SetMinimum(0.01);
      efficiency[ii][0+tightness]->SetMaximum(1.);
      efficiency[ii][0+tightness]->SetMarkerStyle(20);
      efficiency[ii][0+tightness]->SetMarkerSize(1.05);
    }

    // HLT study
    TCanvas c1;
    c1.SetLogy();
    efficiency[0][0+tightness]->SetTitle("");
    efficiency[0][0+tightness]->GetXaxis()->SetTitle("electron #eta");
    efficiency[0][0+tightness]->GetYaxis()->SetTitle("Efficiency");
    efficiency[0][0+tightness]->Draw("PEhist");
    for (int ii=0; ii<3; ii++) {
      efficiency[ii][0+tightness]->Draw("PEhistsame");
    }
    
    TLegend* leg1;
    if (tightness==0 || tightness==1) leg1 = new TLegend(0.15,0.25,0.40,0.4);
    if (tightness==2 || tightness==3) leg1 = new TLegend(0.15,0.75,0.40,0.9);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.03);
    leg1->SetFillColor(0);
    leg1->AddEntry(efficiency[1][0+tightness],"Ele8_CaloIdL_CaloIsoVL");
    leg1->AddEntry(efficiency[2][0+tightness],"Ele17_CaloIdL_CaloIsoVL");
    leg1->AddEntry(efficiency[0][0+tightness],"OR");    
    leg1->Draw();
    
    char tightLevel[10];
    sprintf(tightLevel,"%d",tightness);
    c1.SaveAs("elefake-eta-tightness"+TString(tightLevel)+"_HLTbias.root");
    c1.SaveAs("elefake-eta-tightness"+TString(tightLevel)+"_HLTbias.png");


    // DATA vs MC
    TCanvas c2;
    c2.SetLogy();
    efficiency[0][0+tightness]->SetTitle("");
    efficiency[0][0+tightness]->GetXaxis()->SetTitle("electron #eta");
    efficiency[0][0+tightness]->GetYaxis()->SetTitle("Efficiency");
    efficiency[0][0+tightness]->Draw("PEhist");
    for (int ii=3; ii<6; ii++) {
      efficiency[ii][0+tightness]->Draw("PEhistsame");
    }

    TLegend* leg2;
    if (tightness==0 || tightness==1) leg2 = new TLegend(0.15,0.25,0.40,0.4);
    if (tightness==2 || tightness==3) leg2 = new TLegend(0.15,0.75,0.40,0.9);
    leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.03);
    leg2->SetFillColor(0);
    leg2->AddEntry(efficiency[0][0+tightness],"data");
    leg2->AddEntry(efficiency[3][0+tightness],"QCD");
    leg2->AddEntry(efficiency[4][0+tightness],"Wjets");
    leg2->AddEntry(efficiency[5][0+tightness],"DY");
    leg2->Draw();
    
    sprintf(tightLevel,"%d",tightness);
    c2.SaveAs("elefake-eta-tightness"+TString(tightLevel)+"_datiMc.root");
    c2.SaveAs("elefake-eta-tightness"+TString(tightLevel)+"_datiMc.png");
  }
}


void drawPt() {
  
  TFile *efficiencyFilePt[6];
  // to compare HLT paths
  efficiencyFilePt[0] = TFile::Open("results_data/allData_fakeTeP_may1/Data7TeV/DoubleElectron/4/EleMisidPt.root");
  efficiencyFilePt[1] = TFile::Open("results_data/allData_onlyEle8HLT_fakeTeP_may1/Data7TeV/DoubleElectron/4/EleMisidPt.root");
  efficiencyFilePt[2] = TFile::Open("results_data/allData_onlyEle17HLT_fakeTeP_may1/Data7TeV/DoubleElectron/4/EleMisidPt.root");
  // to compare data with mc                                                                     
  efficiencyFilePt[3] = TFile::Open("macro/fakeTagAndProbe/lastMC_weightedHistos/mergedResultPtQCD.root");
  efficiencyFilePt[4] = TFile::Open("macro/fakeTagAndProbe/lastMC_weightedHistos/mergedResultPtWjets.root");
  efficiencyFilePt[5] = TFile::Open("macro/fakeTagAndProbe/lastMC_weightedHistos/mergedResultPtDY.root");
  
  TH1F *RecoPt[6];
  for (int ii=0; ii<6; ii++) RecoPt[ii] = (TH1F*)efficiencyFilePt[ii]->Get("FakeableJetsPt");

  std::vector<TH1F*> LHIdPt[6];
  for (int ii=0; ii<6; ii++) {
    for (int i=0;i<EgammaLHBasedIDWPs.size();++i) {
      TH1F* aHisto = (TH1F*)efficiencyFilePt[ii]->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt");
      LHIdPt[ii].push_back(aHisto);
    }
  }
  
  EfficiencyEvaluator *ElectronEfficiencyPt[6];
  ElectronEfficiencyPt[0] = new EfficiencyEvaluator("elefake-pt0.root");
  ElectronEfficiencyPt[1] = new EfficiencyEvaluator("elefake-pt1.root");
  ElectronEfficiencyPt[2] = new EfficiencyEvaluator("elefake-pt2.root");
  ElectronEfficiencyPt[3] = new EfficiencyEvaluator("elefake-pt3.root");
  ElectronEfficiencyPt[4] = new EfficiencyEvaluator("elefake-pt4.root");
  ElectronEfficiencyPt[5] = new EfficiencyEvaluator("elefake-pt5.root");

  for (int ii=0; ii<6; ii++) {
    for (int i=0;i<LHIdPt[ii].size();++i) ElectronEfficiencyPt[ii]->AddNumerator(LHIdPt[ii][i]);
    ElectronEfficiencyPt[ii]->SetDenominator(RecoPt[ii]);
    ElectronEfficiencyPt[ii]->ComputeEfficiencies();
  }

  vector<TH1F*> efficiency[6];
  for (int ii=0; ii<6; ii++) efficiency[ii] = ElectronEfficiencyPt[ii]->GetCumulativeEfficiencies();

  for (int tightness=0;tightness<LHIdPt[0].size();++tightness) {

    for (int ii=0; ii<6; ii++) {

      int theColor = ii+1;
      if (ii<3) theColor = ii+1;
      if (ii>2) theColor = ii-1;

      efficiency[ii][0+tightness]->SetLineColor(theColor);
      efficiency[ii][0+tightness]->SetMarkerColor(theColor);
      efficiency[ii][0+tightness]->SetMinimum(0.01);
      efficiency[ii][0+tightness]->SetMaximum(1.);
      efficiency[ii][0+tightness]->SetMarkerStyle(20);
      efficiency[ii][0+tightness]->SetMarkerSize(1.05);
    }

    // HLT study
    TCanvas c1;
    c1.SetLogy();
    efficiency[0][0+tightness]->SetTitle("");
    efficiency[0][0+tightness]->GetXaxis()->SetTitle("electron #pT");
    efficiency[0][0+tightness]->GetYaxis()->SetTitle("Efficiency");
    efficiency[0][0+tightness]->Draw("PEhist");
    for (int ii=0; ii<3; ii++) {
      efficiency[ii][0+tightness]->Draw("PEhistsame");
    }

    TLegend* leg1;
    if (tightness==0 || tightness==1) leg1 = new TLegend(0.15,0.25,0.40,0.4);
    if (tightness==2 || tightness==3) leg1 = new TLegend(0.15,0.75,0.40,0.9);
    leg1->SetFillStyle(0); leg1->SetBorderSize(0); leg1->SetTextSize(0.03);
    leg1->SetFillColor(0);
    leg1->AddEntry(efficiency[1][0+tightness],"Ele8_CaloIdL_CaloIsoVL");
    leg1->AddEntry(efficiency[2][0+tightness],"Ele17_CaloIdL_CaloIsoVL");
    leg1->AddEntry(efficiency[0][0+tightness],"OR");
    leg1->Draw();
    
    char tightLevel[10];
    sprintf(tightLevel,"%d",tightness);
    c1.SaveAs("elefake-pt-tightness"+TString(tightLevel)+"_HLTbias.root");
    c1.SaveAs("elefake-pt-tightness"+TString(tightLevel)+"_HLTbias.png");


    // DATA vs MC
    TCanvas c2;
    c2.SetLogy();
    efficiency[0][0+tightness]->SetTitle("");
    efficiency[0][0+tightness]->GetXaxis()->SetTitle("electron #pT");
    efficiency[0][0+tightness]->GetYaxis()->SetTitle("Efficiency");
    efficiency[0][0+tightness]->Draw("PEhist");
    for (int ii=3; ii<6; ii++) {
      efficiency[ii][0+tightness]->Draw("PEhistsame");
    }

    TLegend* leg2;
    if (tightness==0 || tightness==1) leg2 = new TLegend(0.15,0.25,0.40,0.4);
    if (tightness==2 || tightness==3) leg2 = new TLegend(0.15,0.75,0.40,0.9);
    leg2->SetFillStyle(0); leg2->SetBorderSize(0); leg2->SetTextSize(0.03);
    leg2->SetFillColor(0);
    leg2->AddEntry(efficiency[0][0+tightness],"data");
    leg2->AddEntry(efficiency[3][0+tightness],"QCD");
    leg2->AddEntry(efficiency[4][0+tightness],"Wjets");
    leg2->AddEntry(efficiency[5][0+tightness],"DY");
    leg2->Draw();
    
    sprintf(tightLevel,"%d",tightness);
    c2.SaveAs("elefake-pt-tightness"+TString(tightLevel)+"_datiMc.root");
    c2.SaveAs("elefake-pt-tightness"+TString(tightLevel)+"_datiMc.png");
  }
}
