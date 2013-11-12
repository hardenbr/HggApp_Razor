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

void makePlots(TH1F *signalHistos[2][2][2], TH1F *backgroundHistos[2][2], const char *namevar, const char *axistitle);

int main(int argc, char* argv[]) {

  char inputFileNameSignal[150];
  char inputFileNameBackground[150];
  if ( argc < 3 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: MakeNotePdfPlots fileSignal.root fileBackground.root" << std::endl;
    return 1;
  }
  strcpy(inputFileNameSignal,argv[1]);
  strcpy(inputFileNameBackground,argv[2]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *signalFile = TFile::Open(inputFileNameSignal);
  signalFile->cd("pdfsProducer");

  //[iecal:0=EB;1=EE][iptbin:0=<15GeV;1=>15GeV][iclass:0=nonshowering;1=showering]
  TH1F *dPhiVtxClassEle[2][2][2];
  TH1F *dEtaClassEle[2][2][2];
  TH1F *EoPoutClassEle[2][2][2];
  TH1F *HoEClassEle[2][2][2];
  TH1F *sigmaEtaEtaClassEle[2][2][2];
  TH1F *s9s25ClassEle[2][2][2];

  int iptbin=1; // draw only >15 GeV

  char histo[200];

  for(int iecal=0; iecal<2; iecal++) {
    
    for(int iclass=0; iclass<2; iclass++) {
      
      sprintf(histo,"dPhiVtxClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      dPhiVtxClassEle[iecal][iptbin][iclass]     = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"dEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      dEtaClassEle[iecal][iptbin][iclass]        = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"EoPoutClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      EoPoutClassEle[iecal][iptbin][iclass]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"HoEClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      HoEClassEle[iecal][iptbin][iclass]         = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaEtaEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
      sigmaEtaEtaClassEle[iecal][iptbin][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"s9s25Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
      s9s25ClassEle[iecal][iptbin][iclass] = (TH1F*)gDirectory->Get(histo);

    }

  }


  TFile *backgroundFile = TFile::Open(inputFileNameBackground);
  backgroundFile->cd("pdfsProducer");

  //[iecal:0=EB;1=EE][iptbin:0=<15GeV;1=>15GeV]
  TH1F *dPhiVtxUnsplitJet[2][2];
  TH1F *dEtaUnsplitJet[2][2];
  TH1F *EoPoutUnsplitJet[2][2];
  TH1F *HoEUnsplitJet[2][2];
  TH1F *sigmaEtaEtaUnsplitJet[2][2];
  TH1F *s9s25UnsplitJet[2][2];

  for(int iecal=0; iecal<2; iecal++) {
    
      sprintf(histo,"dPhiVtxUnsplit_electrons_%d_%d",iecal,iptbin);
      dPhiVtxUnsplitJet[iecal][iptbin]     = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"dEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      dEtaUnsplitJet[iecal][iptbin]        = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"EoPoutUnsplit_electrons_%d_%d",iecal,iptbin);
      EoPoutUnsplitJet[iecal][iptbin]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"HoEUnsplit_electrons_%d_%d",iecal,iptbin);
      HoEUnsplitJet[iecal][iptbin]         = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaEtaEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaEtaEtaUnsplitJet[iecal][iptbin] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"s9s25Unsplit_electrons_%d_%d",iecal,iptbin);
      s9s25UnsplitJet[iecal][iptbin] = (TH1F*)gDirectory->Get(histo);

  }

  makePlots(dPhiVtxClassEle,dPhiVtxUnsplitJet,"dPhiIn","#Delta #phi_{in} (rad)");
  makePlots(dEtaClassEle,dEtaUnsplitJet,"dEtaIn","#Delta #eta_{in} (rad)");
  makePlots(EoPoutClassEle,EoPoutUnsplitJet,"EoPout","E_{seed}/p_{out}");
  makePlots(HoEClassEle,HoEUnsplitJet,"HoE","H/E");
  makePlots(sigmaEtaEtaClassEle,sigmaEtaEtaUnsplitJet,"sigmaEtaEta","#sigma_{#eta#eta} (rad^{2})");
  makePlots(s9s25ClassEle,s9s25UnsplitJet,"s9s25","s_{9}/s_{25}");

}

void makePlots(TH1F *signalHistos[2][2][2], TH1F *backgroundHistos[2][2], const char *namevar, const char *axistitle) {

  int iptbin=1;

  for(int iecal=0; iecal<2; iecal++) {

    TCanvas c1;
    c1.SetLogy();

    if(signalHistos[iecal][iptbin][0]) {
      signalHistos[iecal][iptbin][0]->SetLineWidth(2);
      signalHistos[iecal][iptbin][0]->SetLineColor(kRed+1);
      signalHistos[iecal][iptbin][0]->Scale(1.0/signalHistos[iecal][iptbin][0]->Integral());
    }
    if(signalHistos[iecal][iptbin][1]) {
      signalHistos[iecal][iptbin][1]->SetLineWidth(2);
      signalHistos[iecal][iptbin][1]->SetLineColor(kBlue+3);
      signalHistos[iecal][iptbin][1]->Scale(1.0/signalHistos[iecal][iptbin][1]->Integral());
    }
    if(backgroundHistos[iecal][iptbin]) { 
      backgroundHistos[iecal][iptbin]->SetLineWidth(2);
      backgroundHistos[iecal][iptbin]->SetLineColor(kBlack);
      backgroundHistos[iecal][iptbin]->SetLineStyle(kDashed);
      backgroundHistos[iecal][iptbin]->Scale(1.0/backgroundHistos[iecal][iptbin]->Integral());
    }

//     double  max=TMath::Max(signalHistos[iecal][iptbin][0]->GetMaximum(),signalHistos[iecal][iptbin][1]->GetMaximum());
//     max=TMath::Max(max,backgroundHistos[iecal][iptbin]->GetMaximum());
//     signalHistos[iecal][iptbin][0]->SetMaximum(max+sqrt(max));
//     signalHistos[iecal][iptbin][1]->SetMaximum(max+sqrt(max));
//     backgroundHistos[iecal][iptbin]->SetMaximum(max+sqrt(max));

    signalHistos[iecal][iptbin][0]->SetTitle("");
    signalHistos[iecal][iptbin][0]->GetXaxis()->SetTitle(axistitle);
    signalHistos[iecal][iptbin][0]->GetYaxis()->SetTitle("a.u.");
    signalHistos[iecal][iptbin][0]->Draw();
    signalHistos[iecal][iptbin][1]->Draw("same");
    backgroundHistos[iecal][iptbin]->Draw("same");

    TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->AddEntry(signalHistos[iecal][iptbin][0],"non-showering");
    leg->AddEntry(signalHistos[iecal][iptbin][1],"showering");
    leg->AddEntry(backgroundHistos[iecal][iptbin],"jets (all classes)");
    leg->Draw();

    char epsfilename[200];
    char rootfilename[200];
    if(iecal==0) {
      sprintf(epsfilename,"%s_EB.eps",namevar);
      sprintf(rootfilename,"%s_EB.root",namevar);
    }
    else if(iecal==1) {
      sprintf(epsfilename,"%s_EE.eps",namevar); 
      sprintf(rootfilename,"%s_EE.root",namevar);
    }
    c1.SaveAs(epsfilename);
    c1.SaveAs(rootfilename);

  }

}
