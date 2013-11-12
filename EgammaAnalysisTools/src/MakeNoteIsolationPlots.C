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

void makePlots(TH1F *signalHistos[2], TH1F *qcdHistos[2], TH1F* ttbarHistos[2], const char *namevar, const char *axistitle);

int main(int argc, char* argv[]) {

  char inputFileNameSignal[150];
  char inputFileNameQCD[150];
  char inputFileNameTTbar[150];
  if ( argc < 4 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: MakeNotePdfPlots fileSignal.root fileQCD.root fileTTbar" << std::endl;
    return 1;
  }
  strcpy(inputFileNameSignal,argv[1]);
  strcpy(inputFileNameQCD,argv[2]);
  strcpy(inputFileNameTTbar,argv[3]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *signalFile = TFile::Open(inputFileNameSignal);

  // [iecal:0=EB;1=EE]
  TH1F *tkSumPtRelSignal[2];
  TH1F *ecalSumEtRelSignal[2];
  TH1F *hcalSumEtRelSignal[2];
  
  char histo[200];

  for(int iecal=0; iecal<2; iecal++) {

    sprintf(histo,"tkSumPtRel_%d",iecal);
    tkSumPtRelSignal[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"ecalSumEtRel_%d",iecal);
    ecalSumEtRelSignal[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"hcalSumEtRel_%d",iecal);
    hcalSumEtRelSignal[iecal]   = (TH1F*)gDirectory->Get(histo);

  }

  TFile *QCDFile = TFile::Open(inputFileNameQCD);

  // [iecal:0=EB;1=EE]
  TH1F *tkSumPtRelQCD[2];
  TH1F *ecalSumEtRelQCD[2];
  TH1F *hcalSumEtRelQCD[2];
  
  for(int iecal=0; iecal<2; iecal++) {

    sprintf(histo,"tkSumPtRel_%d",iecal);
    tkSumPtRelQCD[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"ecalSumEtRel_%d",iecal);
    ecalSumEtRelQCD[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"hcalSumEtRel_%d",iecal);
    hcalSumEtRelQCD[iecal]   = (TH1F*)gDirectory->Get(histo);

  }

  TFile *TTbarFile = TFile::Open(inputFileNameTTbar);

  // [iecal:0=EB;1=EE]
  TH1F *tkSumPtRelTTbar[2];
  TH1F *ecalSumEtRelTTbar[2];
  TH1F *hcalSumEtRelTTbar[2];
  
  for(int iecal=0; iecal<2; iecal++) {

    sprintf(histo,"tkSumPtRel_%d",iecal);
    tkSumPtRelTTbar[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"ecalSumEtRel_%d",iecal);
    ecalSumEtRelTTbar[iecal]   = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"hcalSumEtRel_%d",iecal);
    hcalSumEtRelTTbar[iecal]   = (TH1F*)gDirectory->Get(histo);

  }

  makePlots(tkSumPtRelSignal,tkSumPtRelQCD,tkSumPtRelTTbar,"tkSumPtRel","#sum p_{T}^{tracks} / p_{T}^{ele}");
  makePlots(ecalSumEtRelSignal,ecalSumEtRelQCD,ecalSumEtRelTTbar,"ecalSumPtRel","#sum E_{T}^{ECAL towers} / E_{T}^{ele}");
  makePlots(hcalSumEtRelSignal,hcalSumEtRelQCD,hcalSumEtRelTTbar,"hcalSumPtRel","#sum E_{T}^{HCAL towers} / E_{T}^{ele}");

}

void makePlots(TH1F *signalHistos[2], TH1F *qcdHistos[2], TH1F *ttbarHistos[2], const char *namevar, const char *axistitle) {

  TCanvas c1;
  c1.SetLogy();

  if(signalHistos[0] && signalHistos[1]) {
    signalHistos[0]->Add(signalHistos[1]);
    signalHistos[0]->SetLineWidth(2);
    signalHistos[0]->SetLineColor(kRed+1);
    signalHistos[0]->Scale(1.0/signalHistos[0]->Integral());
  }
  if(qcdHistos[0] && qcdHistos[1]) {
    qcdHistos[0]->Add(qcdHistos[1]);
    qcdHistos[0]->SetLineWidth(2);
    qcdHistos[0]->SetLineColor(kBlue+1);
    qcdHistos[0]->Scale(1.0/qcdHistos[0]->Integral());
  }
  if(ttbarHistos[0] && ttbarHistos[1]) {
    ttbarHistos[0]->Add(ttbarHistos[1]);
    ttbarHistos[0]->SetLineWidth(2);
    ttbarHistos[0]->SetLineColor(kGreen+1);
    ttbarHistos[0]->Scale(1.0/ttbarHistos[0]->Integral());
  }


  signalHistos[0]->SetTitle("");
  signalHistos[0]->GetXaxis()->SetTitle(axistitle);
  signalHistos[0]->GetYaxis()->SetTitle("a.u.");
  signalHistos[0]->Draw();
  qcdHistos[0]->Draw("same");
  ttbarHistos[0]->Draw("same hist");

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
  leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  leg->AddEntry(signalHistos[0],"signal (Z+jets T & P)");
  leg->AddEntry(qcdHistos[0],"QCD (e.m. enriched + b,c #rightarrow e)");
  leg->AddEntry(ttbarHistos[0],"t #bar{t}");
  leg->Draw();

  char epsfilename[200];
  char rootfilename[200];
  sprintf(epsfilename,"%s_SigBkg.eps",namevar);
  sprintf(rootfilename,"%s_SigBkg.root",namevar);

  c1.SaveAs(epsfilename);
  c1.SaveAs(rootfilename);

}
