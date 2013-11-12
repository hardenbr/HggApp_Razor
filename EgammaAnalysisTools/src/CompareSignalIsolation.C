// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <math.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>

using namespace std;

void makePlots(TH1F *signalHistos1[2], TH1F *signalHistos2[2], const char *namevar, const char *axistitle);

int main(int argc, char* argv[]) {

  char inputFileNameSignal1[150];
  char inputFileNameSignal2[150];
  if ( argc < 3 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: MakeNotePdfPlots fileSignal1.root fileSignal2.root" << std::endl;
    return 1;
  }
  strcpy(inputFileNameSignal1,argv[1]);
  strcpy(inputFileNameSignal2,argv[2]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *signalFile1 = TFile::Open(inputFileNameSignal1);

  // [iecal:0=EB;1=EE]
  TH1F *tkSumPtRel1[2];
  TH1F *ecalSumEtRel1[2];
  TH1F *hcalSumEtRel1[2];
  
  char histo[200];

  for(int iecal=0; iecal<2; iecal++) {

    sprintf(histo,"tkSumPtRel_%d",iecal);
    tkSumPtRel1[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"ecalSumEtRel_%d",iecal);
    ecalSumEtRel1[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"hcalSumEtRel_%d",iecal);
    hcalSumEtRel1[iecal]   = (TH1F*)gDirectory->Get(histo);

  }

  TFile *signalFile2 = TFile::Open(inputFileNameSignal2);

  // [iecal:0=EB;1=EE]
  TH1F *tkSumPtRel2[2];
  TH1F *ecalSumEtRel2[2];
  TH1F *hcalSumEtRel2[2];
  
  for(int iecal=0; iecal<2; iecal++) {

    sprintf(histo,"tkSumPtRel_%d",iecal);
    tkSumPtRel2[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"ecalSumEtRel_%d",iecal);
    ecalSumEtRel2[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"hcalSumEtRel_%d",iecal);
    hcalSumEtRel2[iecal]   = (TH1F*)gDirectory->Get(histo);

  }

  makePlots(tkSumPtRel1,tkSumPtRel2,"tkSumPtRel","#sum p_{T}^{tracks} / p_{T}^{ele}");
  makePlots(ecalSumEtRel1,ecalSumEtRel2,"ecalSumPtRel","#sum E_{T}^{ECAL towers} / E_{T}^{ele}");
  makePlots(hcalSumEtRel1,hcalSumEtRel2,"hcalSumPtRel","#sum E_{T}^{HCAL towers} / E_{T}^{ele}");

}

void makePlots(TH1F *signalHistos1[2], TH1F *signalHistos2[2], const char *namevar, const char *axistitle) {

  TCanvas c1;
  c1.SetLogy();

  for(int iecal=0; iecal<2; iecal++) {

    if(signalHistos1[iecal]) {
      signalHistos1[iecal]->SetLineWidth(2);
      signalHistos1[iecal]->Scale(1.0/signalHistos1[iecal]->Integral());
    }
    if(signalHistos2[iecal]) {
      signalHistos2[iecal]->SetLineWidth(2);
      signalHistos2[iecal]->SetLineStyle(kDashed);
      signalHistos2[iecal]->Scale(1.0/signalHistos2[iecal]->Integral());
    }

  }
    
  signalHistos1[0]->SetLineColor(kRed+1);
  signalHistos2[0]->SetLineColor(kRed+1);
  signalHistos1[1]->SetLineColor(kBlue+1);
  signalHistos2[1]->SetLineColor(kBlue+1);

  signalHistos1[0]->SetTitle("");
  signalHistos1[0]->GetXaxis()->SetTitle(axistitle);
  signalHistos1[0]->GetYaxis()->SetTitle("a.u.");
  signalHistos1[0]->Draw();
  signalHistos1[1]->Draw("same");
  signalHistos2[0]->Draw("same hist");
  signalHistos2[1]->Draw("same hist");

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(signalHistos1[0],"Summer08 barrel");
  leg->AddEntry(signalHistos1[1],"Summer08 endcap");
  leg->AddEntry(signalHistos2[0],"CSA07 barrel");
  leg->AddEntry(signalHistos2[1],"CSA07 endcap");
  leg->Draw();

  char epsfilename[200];
  char rootfilename[200];
  sprintf(epsfilename,"%s.eps",namevar);
  sprintf(rootfilename,"%s.root",namevar);

  c1.SaveAs(epsfilename);
  c1.SaveAs(rootfilename);

}
