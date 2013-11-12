// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TCanvas.h>

using namespace std;

int main(int argc, char* argv[]) {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);


  //************* EB HighPt efficiencies *************
  float eff_EB_CutId_HighPt[5] = { 0.96696 , 0.914958 , 0.852869 , 0.81796 , 0.77299 };
  float eff_EB_LHId_HighPt[5] = { 0.958496 , 0.92175 , 0.865723 , 0.816584 , 0.774055 };
  //************* EE HighPt efficiencies *************
  float eff_EE_CutId_HighPt[5] = { 0.936028 , 0.826423 , 0.749348 , 0.678982 , 0.598987 };
  float eff_EE_LHId_HighPt[5] = { 0.935676 , 0.879968 , 0.763354 , 0.685813 , 0.608444 };
  //************* EB LowPt efficiencies *************
  float eff_EB_CutId_LowPt[5] = { 0.871968 , 0.735454 , 0.650924 , 0.608109 , 0.551101 };
  float eff_EB_LHId_LowPt[5] = { 0.880872 , 0.807528 , 0.725145 , 0.650965 , 0.583139 };
  //************* EE LowPt efficiencies *************
  float eff_EE_CutId_LowPt[5] = { 0.71685 , 0.58594 , 0.469225 , 0.394948 , 0.336686 };
  float eff_EE_LHId_LowPt[5] = { 0.819503 , 0.700758 , 0.539236 , 0.413665 , 0.298364 };
  //************* EB HighPt fakes *************
  float fr_EB_CutId_HighPt[5] = { 0.111561 , 0.0841345 , 0.0571164 , 0.0497258 , 0.0441188 };
  float fr_EB_LHId_HighPt[5] = { 0.0866837 , 0.0642953 , 0.0505614 , 0.04346 , 0.0394024 };
  //************* EE HighPt fakes *************
  float fr_EE_CutId_HighPt[5] = { 0.155972 , 0.0964785 , 0.0655005 , 0.0459625 , 0.0376067 };
  float fr_EE_LHId_HighPt[5] = { 0.151757 , 0.109348 , 0.0717574 , 0.0481252 , 0.0364888 };
  //************* EB LowPt fakes *************
  float fr_EB_CutId_LowPt[5] = { 0.0740105 , 0.0459741 , 0.0351585 , 0.0304999 , 0.0271889 };
  float fr_EB_LHId_LowPt[5] = { 0.0683885 , 0.0485206 , 0.0371129 , 0.0301147 , 0.0261233 };
  //************* EE LowPt fakes *************
  float fr_EE_CutId_LowPt[5] = { 0.118747 , 0.0784185 , 0.0492684 , 0.0330966 , 0.026463 };
  float fr_EE_LHId_LowPt[5] = { 0.153119 , 0.101379 , 0.0579594 , 0.0314032 , 0.0194664 };

/* OLD ROC with no OneOverEMinusOneOverP
// MC efficiencies on W+jets Spring 11 (EB)
float eff_EB_CutId_HighPt[5] = { 0.97, 0.91, 0.85, 0.82, 0.77 };
float eff_EB_LHId_HighPt[5] = { 0.96, 0.93, 0.87, 0.80, 0.73 };

float eff_EB_CutId_LowPt[5] = { 0.87, 0.74, 0.65, 0.61, 0.55 };
float eff_EB_LHId_LowPt[5] = { 0.90, 0.83, 0.73, 0.62, 0.50 };

// MC efficiencies on W+jets Spring 11 (EE)
float eff_EE_CutId_HighPt[5] = { 0.94, 0.83, 0.75, 0.68, 0.60 };
float eff_EE_LHId_HighPt[5] = { 0.95, 0.90, 0.78, 0.68, 0.54 };

float eff_EE_CutId_LowPt[5] = { 0.72, 0.59, 0.47, 0.39, 0.34 };
float eff_EE_LHId_LowPt[5] = { 0.84, 0.73, 0.54, 0.37, 0.21 };

// MC fake rates on W+jets Spring 11 (EB)
float fr_EB_CutId_HighPt[5] = { 0.112, 0.084, 0.057, 0.050, 0.044 };
float fr_EB_LHId_HighPt[5] = { 0.102, 0.075, 0.056, 0.046, 0.039 };

float fr_EB_CutId_LowPt[5] = { 0.074, 0.046, 0.035, 0.030, 0.027 };
float fr_EB_LHId_LowPt[5] = { 0.088, 0.055, 0.040, 0.030, 0.022 };

// MC fake rates on W+jets Spring 11 (EE)
float fr_EE_CutId_HighPt[5] = { 0.156, 0.096, 0.066, 0.046, 0.038 };
float fr_EE_LHId_HighPt[5] = { 0.163, 0.119, 0.075, 0.047, 0.031 };

float fr_EE_CutId_LowPt[5] = { 0.119, 0.078, 0.049, 0.033, 0.026 };
float fr_EE_LHId_LowPt[5] = { 0.171, 0.109, 0.062, 0.033, 0.015 };
*/
  

  TGraph *ROC_EB_CutId_HighPt = new TGraph(5, eff_EB_CutId_HighPt, fr_EB_CutId_HighPt);
  TGraph *ROC_EB_LHId_HighPt = new TGraph(5, eff_EB_LHId_HighPt, fr_EB_LHId_HighPt);

  TGraph *ROC_EE_CutId_HighPt = new TGraph(5, eff_EE_CutId_HighPt, fr_EE_CutId_HighPt);
  TGraph *ROC_EE_LHId_HighPt = new TGraph(5, eff_EE_LHId_HighPt, fr_EE_LHId_HighPt);

  TGraph *ROC_EB_CutId_LowPt = new TGraph(5, eff_EB_CutId_LowPt, fr_EB_CutId_LowPt);
  TGraph *ROC_EB_LHId_LowPt = new TGraph(5, eff_EB_LHId_LowPt, fr_EB_LHId_LowPt);

  TGraph *ROC_EE_CutId_LowPt = new TGraph(5, eff_EE_CutId_LowPt, fr_EE_CutId_LowPt);
  TGraph *ROC_EE_LHId_LowPt = new TGraph(5, eff_EE_LHId_LowPt, fr_EE_LHId_LowPt);

  // high pT ROC EB
  TCanvas *highPtEB = new TCanvas("highPtEB","highPtEB",600,400);
  ROC_EB_CutId_HighPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_EB_CutId_HighPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_EB_CutId_HighPt->SetMarkerStyle(20);
  ROC_EB_CutId_HighPt->SetMarkerSize(1.05);
  ROC_EB_LHId_HighPt->SetMarkerStyle(21);
  ROC_EB_LHId_HighPt->SetMarkerSize(1.05);
  ROC_EB_CutId_HighPt->SetLineColor(1);
  ROC_EB_CutId_HighPt->SetMarkerColor(1);
  ROC_EB_LHId_HighPt->SetLineColor(2);
  ROC_EB_LHId_HighPt->SetMarkerColor(2);

  TLegend* leg = new TLegend(0.35,0.25,0.60,0.50);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(ROC_EB_CutId_HighPt, "Cut Id","pl"); 
  leg->AddEntry(ROC_EB_LHId_HighPt, "LH Id","pl"); 

  ROC_EB_CutId_HighPt->Draw("apl");
  ROC_EB_LHId_HighPt->Draw("pl");
  leg->Draw();
  highPtEB->SaveAs("roc_EB_HighPt.eps");
  highPtEB->SaveAs("roc_EB_HighPt.root");
  highPtEB->SaveAs("roc_EB_HighPt.png");

  // high pT ROC EE
  TCanvas *highPtEE = new TCanvas("highPtEE","highPtEE",600,400);
  ROC_EE_CutId_HighPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_EE_CutId_HighPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_EE_CutId_HighPt->SetMarkerStyle(20);
  ROC_EE_CutId_HighPt->SetMarkerSize(1.05);
  ROC_EE_LHId_HighPt->SetMarkerStyle(21);
  ROC_EE_LHId_HighPt->SetMarkerSize(1.05);
  ROC_EE_CutId_HighPt->SetLineColor(1);
  ROC_EE_CutId_HighPt->SetMarkerColor(1);
  ROC_EE_LHId_HighPt->SetLineColor(2);
  ROC_EE_LHId_HighPt->SetMarkerColor(2);

  ROC_EE_CutId_HighPt->Draw("apl");
  ROC_EE_LHId_HighPt->Draw("pl");
  leg->Draw();
  highPtEE->SaveAs("roc_EE_HighPt.eps");
  highPtEE->SaveAs("roc_EE_HighPt.root");
  highPtEE->SaveAs("roc_EE_HighPt.png");

  // low pT ROC EB
  TCanvas *lowPtEB = new TCanvas("lowPtEB","lowPtEB",600,400);
  ROC_EB_CutId_LowPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_EB_CutId_LowPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_EB_CutId_LowPt->SetMarkerStyle(20);
  ROC_EB_CutId_LowPt->SetMarkerSize(1.05);
  ROC_EB_LHId_LowPt->SetMarkerStyle(21);
  ROC_EB_LHId_LowPt->SetMarkerSize(1.05);
  ROC_EB_CutId_LowPt->SetLineColor(1);
  ROC_EB_CutId_LowPt->SetMarkerColor(1);
  ROC_EB_LHId_LowPt->SetLineColor(2);
  ROC_EB_LHId_LowPt->SetMarkerColor(2);

  ROC_EB_CutId_LowPt->Draw("apl");
  ROC_EB_LHId_LowPt->Draw("pl");
  leg->Draw();
  lowPtEB->SaveAs("roc_EB_LowPt.eps");
  lowPtEB->SaveAs("roc_EB_LowPt.root");
  lowPtEB->SaveAs("roc_EB_LowPt.png");

  // low pT ROC EE
  TCanvas *lowPtEE = new TCanvas("lowPtEE","lowPtEE",600,400);
  ROC_EE_CutId_LowPt->GetXaxis()->SetTitle("Signal efficiency");
  ROC_EE_CutId_LowPt->GetYaxis()->SetTitle("Background efficiency");
  ROC_EE_CutId_LowPt->SetMarkerStyle(20);
  ROC_EE_CutId_LowPt->SetMarkerSize(1.05);
  ROC_EE_LHId_LowPt->SetMarkerStyle(21);
  ROC_EE_LHId_LowPt->SetMarkerSize(1.05);
  ROC_EE_CutId_LowPt->SetLineColor(1);
  ROC_EE_CutId_LowPt->SetMarkerColor(1);
  ROC_EE_LHId_LowPt->SetLineColor(2);
  ROC_EE_LHId_LowPt->SetMarkerColor(2);

  ROC_EE_CutId_LowPt->Draw("apl");
  ROC_EE_LHId_LowPt->Draw("pl");
  leg->Draw();
  lowPtEE->SaveAs("roc_EE_LowPt.eps");
  lowPtEE->SaveAs("roc_EE_LowPt.root");
  lowPtEE->SaveAs("roc_EE_LowPt.png");


}
