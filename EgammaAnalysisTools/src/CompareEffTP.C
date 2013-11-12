#include "TCanvas.h"
#include "TFile.h"
#include "RooHist.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TAxis.h"
#include "TString.h"
#include <vector>
#include <sstream>
#include <fstream>

void compare(const char* file1, const char* file2, const char* variable="vertices", const char* effPoint="80", bool writeEBEE=false);
float eff_Cut[2], eff_LH[2];

int doEff() {

  // for ROC curve
  std::ostringstream effCutEB;
  std::ostringstream effCutEE;
  std::ostringstream effLHEB;
  std::ostringstream effLHEE;

  std::vector<TString> CutFile, LHFile, WP;
  CutFile.push_back("WP95"); LHFile.push_back("LHVeryLoose"); WP.push_back("95");
  CutFile.push_back("WP90"); LHFile.push_back("LHLoose"); WP.push_back("90");
  CutFile.push_back("WP85"); LHFile.push_back("LHMedium"); WP.push_back("85");
  CutFile.push_back("WP80"); LHFile.push_back("LHTight"); WP.push_back("80");
  CutFile.push_back("WP70"); LHFile.push_back("LHHyperTight"); WP.push_back("70");

  effCutEB << "float eff_EB_CutId_" << "[ " << WP.size() << " ] = { ";
  effCutEE << "float eff_EE_CutId_" << "[ " << WP.size() << " ] = { ";

  effLHEB << "float eff_EB_LHId_" << "[ " << WP.size() << " ] = { ";
  effLHEE << "float eff_EE_LHId_" << "[ " << WP.size() << " ] = { ";

  // eta
  for(int i=0; i<(int)WP.size(); i++) {
    TString CutFilePath = TString("NewCatOpt_1oEm1oP/HighPt/effData_Eta_")+CutFile[i]+TString(".root");
    TString LHFilePath = TString("NewCatOpt_1oEm1oP/HighPt/effData_Eta_")+LHFile[i]+TString(".root");
    compare(CutFilePath.Data(),LHFilePath.Data(),"eta",WP[i].Data(),true);
    effCutEB << eff_Cut[0] ;
      effCutEE << eff_Cut[1] ;
      effLHEB << eff_LH[0];
      effLHEE << eff_LH[1] ;
    if(i==(int)WP.size()-1) {
      effCutEB << " };"  ;
      effCutEE << " };"  ;
      effLHEB << " };"  ;
      effLHEE << " };"  ;
    } else {
      effCutEB << " , " ;
      effCutEE << " , " ;
      effLHEB << " , " ;
      effLHEE << " , " ;
    }
  }

  std::ostringstream fileName;
  fileName << "efficiencies.txt";

  std::ofstream out;
  out.open(fileName.str().c_str());
  out << "//************* EB efficiencies *************" << std::endl;
  out << effCutEB.str() << std::endl;
  out << effLHEB.str() << std::endl;

  out << "//************* EE efficiencies *************" << std::endl;
  out << effCutEE.str() << std::endl;
  out << effLHEE.str() << std::endl;
  

  // pT
  for(int i=0; i<(int)WP.size(); i++) {
    TString CutFilePath = TString("NewCatOpt_1oEm1oP/HighPt/effData_Pt_")+CutFile[i]+TString(".root");
    TString LHFilePath = TString("NewCatOpt_1oEm1oP/HighPt/effData_Pt_")+LHFile[i]+TString(".root");
    compare(CutFilePath.Data(),LHFilePath.Data(),"pt",WP[i].Data());
  }
  // vertices
  for(int i=0; i<(int)WP.size(); i++) {
    TString CutFilePath = TString("NewCatOpt_1oEm1oP/HighPt/effData_Vertices_")+CutFile[i]+TString(".root");
    TString LHFilePath = TString("NewCatOpt_1oEm1oP/HighPt/effData_Vertices_")+LHFile[i]+TString(".root");
    compare(CutFilePath.Data(),LHFilePath.Data(),"vertices",WP[i].Data());
  }

  return 0;
}

void compare(const char* file1, const char* file2, const char* variable, const char* effPoint, bool writeEBEE) {

  gROOT->SetStyle("Plain");

  TString dir = TString("eleIDdir/")+TString(variable)+TString("/fit_eff_plots/")+TString(variable)+TString("_PLOT");

  TFile *TFile1 = TFile::Open(file1);
  TCanvas *c1 = (TCanvas*)TFile1->Get(dir.Data());
  RooHist* h1 = (RooHist*)c1->GetPrimitive("hxy_fit_eff");
  
  TFile *TFile2 = TFile::Open(file2);
  TCanvas *c2 = (TCanvas*)TFile2->Get(dir.Data());
  RooHist* h2 = (RooHist*)c2->GetPrimitive("hxy_fit_eff");

  TCanvas *c3 = new TCanvas("c3","c3",600,400);
  h1->GetXaxis()->SetTitle(variable);
  h2->SetMarkerColor(kRed);
  h2->SetLineColor(kRed);
  h1->SetMinimum(0.0);
  h1->SetMaximum(1.0);
  h1->Draw("ap");
  h2->Draw("p");

  // compute the weighted average
  Double_t* x1 = h1->GetX();
  Double_t* y1 = h1->GetY();
  Double_t* y2 = h2->GetY();
  Double_t* ey1 = h1->GetEYlow();
  Double_t* ey2 = h2->GetEYlow();
  Int_t nbins = h1->GetN();

  Double_t num1, denom1, num2, denom2;
  num1 = denom1 = num2 = denom2 = 0.0;

  Double_t num1_eta[2], denom1_eta[2], num2_eta[2], denom2_eta[2];
  num1_eta[0] = denom1_eta[0] = num2_eta[0] = denom2_eta[0] = 0.0;
  num1_eta[1] = denom1_eta[1] = num2_eta[1] = denom2_eta[1] = 0.0;

  for(int i=0; i<nbins; i++) {
    std::cout << "i = " << i << " WP eff = " << ey1[i] << " LH eff = " << ey2[i] << std::endl;
    if(fabs(ey1[i])>0) {
      num1 += y1[i]/ey1[i]/ey1[i];
      denom1 += 1/ey1[i]/ey1[i];

      if(writeEBEE) {
	// now compute the EB/EE separately
	int subdet;
	if(fabs(x1[i])<1.479) subdet = 0;
	else subdet = 1;
	num1_eta[subdet] += y1[i]/ey1[i]/ey1[i];
	denom1_eta[subdet] += 1/ey1[i]/ey1[i];
      }
    }
    if(fabs(ey2[i])>0) {
      num2 += y2[i]/ey2[i]/ey2[i];
      denom2 += 1/ey2[i]/ey2[i];

      if(writeEBEE) {
	// now compute the EB/EE separately
	int subdet;
	if(fabs(x1[i])<1.479) subdet = 0;
	else subdet = 1;
	num2_eta[subdet] += y2[i]/ey2[i]/ey2[i];
	denom2_eta[subdet] += 1/ey2[i]/ey2[i];
      }
    }
  }
  // std::cout << "EFFICIENCY FOR " << effPoint << " (WP) = " << num1/denom1 
  // 	    << " (LH) = " << num2/denom2 << std::endl;

  char CutIdEff[50], LHIdEff[50];
  char CutIdEff_EB[50], LHIdEff_EB[50], CutIdEff_EE[50], LHIdEff_EE[50];
  
  sprintf(CutIdEff,"Cut Id WP%s (#epsilon = %.2f)", effPoint, num1/denom1);
  sprintf(LHIdEff,"LH Id WP%s (#epsilon = %.2f)", effPoint, num2/denom2);
  
  if(writeEBEE) {
    sprintf(CutIdEff_EB,"Cut Id WP%s (#epsilon EB = %.2f)", effPoint, num1_eta[0]/denom1_eta[0]);
    sprintf(CutIdEff_EE,"Cut Id WP%s (#epsilon EE = %.2f)", effPoint, num1_eta[1]/denom1_eta[1]);

    sprintf(LHIdEff_EB,"LH Id WP%s (#epsilon EB = %.2f)", effPoint, num2_eta[0]/denom2_eta[0]);
    sprintf(LHIdEff_EE,"LH Id WP%s (#epsilon EE = %.2f)", effPoint, num2_eta[1]/denom2_eta[1]);
  }

  TLegend* leg = new TLegend(0.15, 0.20, 0.85, 0.4);
  leg->SetFillStyle(0); leg->SetBorderSize(0.); leg->SetTextSize(0.025);
  leg->SetFillColor(0);
  if(writeEBEE) {
    leg->AddEntry(h1,CutIdEff_EB,"p");
    leg->AddEntry(h1,CutIdEff_EE,"p");
    leg->AddEntry(h2,LHIdEff_EB,"p");
    leg->AddEntry(h2,LHIdEff_EE,"p");
  } else {
    leg->AddEntry(h1,CutIdEff,"p");
    leg->AddEntry(h2,LHIdEff,"p");
  }
  leg->Draw();

  TString fileOut = TString("eff_vs_")+TString(variable)+TString("_")+TString(effPoint)+TString("_HighPt.png");
  c3->SaveAs(fileOut.Data());

  fileOut = TString("eff_vs_")+TString(variable)+TString("_")+TString(effPoint)+TString("_HighPt.eps");
  c3->SaveAs(fileOut.Data());

  if(writeEBEE) {
    for(int iecal=0;iecal<2;iecal++) {
      eff_Cut[iecal] = num1_eta[iecal]/denom1_eta[iecal];
      eff_LH[iecal] = num2_eta[iecal]/denom2_eta[iecal];
    }
  }

  return;
}
