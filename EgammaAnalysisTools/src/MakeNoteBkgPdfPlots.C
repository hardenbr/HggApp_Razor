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

void makePlots(TH1F *qcdHistos[2][2], TH1F *wjetsHistos[2][2], const char *namevar, const char *axistitle);

int main(int argc, char* argv[]) {

  char inputFileNameQCD[150];
  char inputFileNameWjets[150];
  if ( argc < 3 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: MakeNotePdfPlots fileQCD.root fileWjets.root" << std::endl;
    return 1;
  }
  strcpy(inputFileNameQCD,argv[1]);
  strcpy(inputFileNameWjets,argv[2]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *qcdFile = TFile::Open(inputFileNameQCD);
  qcdFile->cd("");

  //[iecal:0=EB;1=EE][iptbin:0=<15GeV;1=>15GeV]
  TH1F *dPhiVtxUnsplitQCD[2][2];
  TH1F *dEtaUnsplitQCD[2][2];
  TH1F *EoPoutUnsplitQCD[2][2];
  TH1F *HoEUnsplitQCD[2][2];
  TH1F *sigmaEtaEtaUnsplitQCD[2][2];
  TH1F *s9s25UnsplitQCD[2][2];

  int iptbin=1; // draw only >15 GeV

  char histo[200];

  for(int iecal=0; iecal<2; iecal++) {
    
    sprintf(histo,"dPhiVtxUnsplit_electrons_%d_%d",iecal,iptbin);
    dPhiVtxUnsplitQCD[iecal][iptbin]     = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"dEtaUnsplit_electrons_%d_%d",iecal,iptbin);
    dEtaUnsplitQCD[iecal][iptbin]        = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"EoPoutUnsplit_electrons_%d_%d",iecal,iptbin);
    EoPoutUnsplitQCD[iecal][iptbin]      = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"HoEUnsplit_electrons_%d_%d",iecal,iptbin);
    HoEUnsplitQCD[iecal][iptbin]         = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"sigmaEtaEtaUnsplit_electrons_%d_%d",iecal,iptbin);
    sigmaEtaEtaUnsplitQCD[iecal][iptbin] = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"s9s25Unsplit_electrons_%d_%d",iecal,iptbin);
    s9s25UnsplitQCD[iecal][iptbin] = (TH1F*)gDirectory->Get(histo);

  }


  TFile *wjetsFile = TFile::Open(inputFileNameWjets);
  wjetsFile->cd("pdfsProducer");

  //[iecal:0=EB;1=EE][iptbin:0=<15GeV;1=>15GeV]
  TH1F *dPhiVtxUnsplitWjets[2][2];
  TH1F *dEtaUnsplitWjets[2][2];
  TH1F *EoPoutUnsplitWjets[2][2];
  TH1F *HoEUnsplitWjets[2][2];
  TH1F *sigmaEtaEtaUnsplitWjets[2][2];
  TH1F *s9s25UnsplitWjets[2][2];

  for(int iecal=0; iecal<2; iecal++) {
    
      sprintf(histo,"dPhiVtxUnsplit_electrons_%d_%d",iecal,iptbin);
      dPhiVtxUnsplitWjets[iecal][iptbin]     = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"dEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      dEtaUnsplitWjets[iecal][iptbin]        = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"EoPoutUnsplit_electrons_%d_%d",iecal,iptbin);
      EoPoutUnsplitWjets[iecal][iptbin]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"HoEUnsplit_electrons_%d_%d",iecal,iptbin);
      HoEUnsplitWjets[iecal][iptbin]         = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaEtaEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaEtaEtaUnsplitWjets[iecal][iptbin] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"s9s25Unsplit_electrons_%d_%d",iecal,iptbin);
      s9s25UnsplitWjets[iecal][iptbin] = (TH1F*)gDirectory->Get(histo);

  }

  makePlots(dPhiVtxUnsplitQCD,dPhiVtxUnsplitWjets,"dPhiIn","#Delta #phi_{in} (rad)");
  makePlots(dEtaUnsplitQCD,dEtaUnsplitWjets,"dEtaIn","#Delta #eta_{in} (rad)");
  makePlots(EoPoutUnsplitQCD,EoPoutUnsplitWjets,"EoPout","E_{seed}/p_{out}");
  makePlots(HoEUnsplitQCD,HoEUnsplitWjets,"HoE","H/E");
  makePlots(sigmaEtaEtaUnsplitQCD,sigmaEtaEtaUnsplitWjets,"sigmaEtaEta","#sigma_{#eta#eta} (rad^{2})");
  makePlots(s9s25UnsplitQCD,s9s25UnsplitWjets,"s9s25","s_{9}/s_{25}");

}

void makePlots(TH1F *qcdHistos[2][2], TH1F *wjetsHistos[2][2], const char *namevar, const char *axistitle) {

  int iptbin=1;

  for(int iecal=0; iecal<2; iecal++) {

    TCanvas c1;
    c1.SetLogy();

    cout << "iecal = " << iecal << " iptbin = " << iptbin 
         << " qcdHistos[iecal][iptbin] = " << qcdHistos[iecal][iptbin] <<  endl;

    if(qcdHistos[iecal][iptbin]) {
      qcdHistos[iecal][iptbin]->Rebin(2);
      qcdHistos[iecal][iptbin]->SetLineWidth(2);
      qcdHistos[iecal][iptbin]->SetLineColor(kRed+1);
      qcdHistos[iecal][iptbin]->Scale(1.0/qcdHistos[iecal][iptbin]->Integral());
    }
    if(wjetsHistos[iecal][iptbin]) { 
      wjetsHistos[iecal][iptbin]->Rebin(2);
      wjetsHistos[iecal][iptbin]->SetLineWidth(2);
      wjetsHistos[iecal][iptbin]->SetLineColor(kBlue+1);
      wjetsHistos[iecal][iptbin]->Scale(1.0/wjetsHistos[iecal][iptbin]->Integral());
    }

//     double  max=TMath::Max(signalHistos[iecal][iptbin][0]->GetMaximum(),signalHistos[iecal][iptbin][1]->GetMaximum());
//     max=TMath::Max(max,backgroundHistos[iecal][iptbin]->GetMaximum());
//     signalHistos[iecal][iptbin][0]->SetMaximum(max+sqrt(max));
//     signalHistos[iecal][iptbin][1]->SetMaximum(max+sqrt(max));
//     backgroundHistos[iecal][iptbin]->SetMaximum(max+sqrt(max));

    qcdHistos[iecal][iptbin]->SetTitle("");
    qcdHistos[iecal][iptbin]->GetXaxis()->SetTitle(axistitle);
    qcdHistos[iecal][iptbin]->GetYaxis()->SetTitle("a.u.");
    qcdHistos[iecal][iptbin]->Draw();
    wjetsHistos[iecal][iptbin]->Draw("same");

    TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->AddEntry(qcdHistos[iecal][iptbin],"jets from QCD");
    leg->AddEntry(wjetsHistos[iecal][iptbin],"jets from W+jets");
    leg->Draw();

    char epsfilename[200];
    char rootfilename[200];
    if(iecal==0) {
      sprintf(epsfilename,"%s_qcdVsWjets_EB.eps",namevar);
      sprintf(rootfilename,"%s_qcdVsWjets_EB.root",namevar);
    }
    else if(iecal==1) {
      sprintf(epsfilename,"%s_qcdVsWjets_EE.eps",namevar); 
      sprintf(rootfilename,"%s_qcdVsWjets_EE.root",namevar);
    }
    c1.SaveAs(epsfilename);
    c1.SaveAs(rootfilename);

  }

}
