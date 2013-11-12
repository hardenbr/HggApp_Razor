// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ROOT includes
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TString.h>
#include <TAxis.h>

using namespace std;
int MODE;

void makePlots(TH1F *signalHistos1[2][2], TH1F *signalHistos2[2][2], const char *namevar, const char *axistitle);

int main(int argc, char* argv[]) {

  char inputFileNameSignal1[150];
  char inputFileNameSignal2[150];
  if ( argc < 4 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: MakeNotePdfPlots fileSignalMC.root fileSignalDATA.root 0" << std::endl;
    std::cout << "third argument = 0 ==> signal, MC vs data. 1 ==> signal vs bkg" << std::endl;
    return 1;
  }
  strcpy(inputFileNameSignal1,argv[1]);
  strcpy(inputFileNameSignal2,argv[2]);
  MODE = atoi(argv[3]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *signalFile1 = TFile::Open(inputFileNameSignal1);
  signalFile1->cd();

  TH1F *dPhiClassEle1[2][2];
  TH1F *dEtaClassEle1[2][2];
  TH1F *EoPClassEle1[2][2];
  TH1F *HoEClassEle1[2][2];
  TH1F *sigmaIEtaIEtaClassEle1[2][2];
  TH1F *sigmaIPhiIPhiClassEle1[2][2];
  TH1F *fbremClassEle1[2][2];
  TH1F *lhClassEle1[2][2];

  char histo[200];

  int iptbin = 1; // > 15 GeV
  for(int iecal=0; iecal<2; iecal++) {

    // iclass = 0: 0 - brem clusters
    // iclass = 1: >=1 - brem clusters
    for(int iclass=0; iclass<2; iclass++) {
      
      sprintf(histo,"dPhiClass_electrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
      dPhiClassEle1[iecal][iclass]     = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"dEtaClass_electrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
      dEtaClassEle1[iecal][iclass]        = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"EoPClass_electrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
      EoPClassEle1[iecal][iclass]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"HoEClass_electrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
      HoEClassEle1[iecal][iclass]         = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaIEtaIEtaClass_electrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
      sigmaIEtaIEtaClassEle1[iecal][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaIPhiIPhiClass_electrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
      sigmaIPhiIPhiClassEle1[iecal][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"fBremClass_electrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
      fbremClassEle1[iecal][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"lhClass_electrons_subdet%d_ptbin%d_class%d",iecal,iptbin,iclass);
      lhClassEle1[iecal][iclass] = (TH1F*)gDirectory->Get(histo); 
    
      dPhiClassEle1[iecal][iclass]   -> SetMinimum(0.001);
      dEtaClassEle1[iecal][iclass]   -> SetMinimum(0.001);
      EoPClassEle1[iecal][iclass]    -> SetMinimum(0.001);
      HoEClassEle1[iecal][iclass]   -> SetMinimum(0.001);
      sigmaIEtaIEtaClassEle1[iecal][iclass] -> SetMinimum(0.1);
      sigmaIPhiIPhiClassEle1[iecal][iclass] -> SetMinimum(0.1);
      fbremClassEle1[iecal][iclass] -> SetMinimum(0.001);
      lhClassEle1[iecal][iclass] -> SetMinimum(0.001);
    
    }
  }

  TFile *signalFile2 = TFile::Open(inputFileNameSignal2);
  signalFile2->cd();

  TH1F *dPhiClassEle2[2][2];
  TH1F *dEtaClassEle2[2][2];
  TH1F *EoPClassEle2[2][2];
  TH1F *HoEClassEle2[2][2];
  TH1F *sigmaIEtaIEtaClassEle2[2][2];
  TH1F *sigmaIPhiIPhiClassEle2[2][2];
  TH1F *fbremClassEle2[2][2];
  TH1F *lhClassEle2[2][2];

  for(int iecal=0; iecal<2; iecal++) {

    // iclass = 0: 0 - brem clusters
    // iclass = 1: >=1 - brem clusters
    for(int iclass=0; iclass<2; iclass++) {
      
      char hypothesis[200];
      if(MODE==0) sprintf(hypothesis,"electrons");
      else sprintf(hypothesis,"hadrons");
      
      sprintf(histo,"dPhiClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
      dPhiClassEle2[iecal][iclass]     = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"dEtaClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
      dEtaClassEle2[iecal][iclass]        = (TH1F*)gDirectory->Get(histo);
      sprintf(histo,"EoPClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
      EoPClassEle2[iecal][iclass]      = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"HoEClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
      HoEClassEle2[iecal][iclass]         = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaIEtaIEtaClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
      sigmaIEtaIEtaClassEle2[iecal][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"sigmaIPhiIPhiClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
      sigmaIPhiIPhiClassEle2[iecal][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"fBremClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
      fbremClassEle2[iecal][iclass] = (TH1F*)gDirectory->Get(histo); 
      sprintf(histo,"lhClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
      lhClassEle2[iecal][iclass] = (TH1F*)gDirectory->Get(histo); 

      dPhiClassEle2[iecal][iclass]   -> SetMinimum(0.001);
      dEtaClassEle2[iecal][iclass]   -> SetMinimum(0.001);
      EoPClassEle2[iecal][iclass]    -> SetMinimum(0.001);
      HoEClassEle2[iecal][iclass]   -> SetMinimum(0.001);
      sigmaIEtaIEtaClassEle2[iecal][iclass] -> SetMinimum(0.1);
      sigmaIPhiIPhiClassEle2[iecal][iclass] -> SetMinimum(0.1);
      fbremClassEle2[iecal][iclass]   -> SetMinimum(0.001);
      lhClassEle2[iecal][iclass] -> SetMinimum(0.001);

    }
  }

  makePlots(dPhiClassEle1,dPhiClassEle2,"dPhiIn","#Delta #phi_{in} (rad)");
  makePlots(dEtaClassEle1,dEtaClassEle2,"dEtaIn","#Delta #eta_{in} (rad)");
  makePlots(EoPClassEle1,EoPClassEle2,"EoP","E_{SC}/p_{in}");
  makePlots(HoEClassEle1,HoEClassEle2,"HoE","H/E");
  makePlots(sigmaIEtaIEtaClassEle1,sigmaIEtaIEtaClassEle2,"sigmaIEtaIEta","#sigma_{#eta#eta} (rad)");
  makePlots(sigmaIPhiIPhiClassEle1,sigmaIPhiIPhiClassEle2,"sigmaIPhiIPhi","#sigma_{#phi#phi} (rad)");
  makePlots(fbremClassEle1,fbremClassEle2,"fbrem","f_{brem}");
  makePlots(lhClassEle1,lhClassEle2,"lh","LH output");

}

void makePlots(TH1F *signalHistos1[2][2], TH1F *signalHistos2[2][2], const char *namevar, const char *axistitle) {

  TString intLumi="2.88";

  int iptbin=1;

  for(int iecal=0; iecal<2; iecal++) {

    TCanvas c1;
    // c1.SetLogy();

    // iclass = 0: 0 - brem clusters
    // iclass = 1: >=1 - brem clusters
    for(int iclass=0; iclass<2; iclass++) {

    // data (MODE=0) or bkg (MODE=1)
    if(signalHistos2[iecal][iclass]) {
      signalHistos2[iecal][iclass]->Rebin(5);
      signalHistos2[iecal][iclass]->SetLineWidth(2);
      if(MODE==0) {
        signalHistos2[iecal][iclass]->SetLineColor(kRed+1+3*iclass);
        signalHistos2[iecal][iclass]->SetMarkerColor(kRed+1+3*iclass);
      } else {
        signalHistos2[iecal][iclass]->SetLineColor(kViolet+1+3*iclass);
        signalHistos2[iecal][iclass]->SetMarkerColor(kViolet+1+3*iclass);
      }
      signalHistos2[iecal][iclass]->SetMarkerStyle(8);
      signalHistos2[iecal][iclass]->SetMarkerSize(1.0);
    }

    // MC
    if(signalHistos1[iecal][iclass]) {
      signalHistos1[iecal][iclass]->Rebin(5);
      signalHistos1[iecal][iclass]->SetLineWidth(2);
      signalHistos1[iecal][iclass]->SetLineColor(kRed+1+3*iclass);
      signalHistos1[iecal][iclass]->Scale((float)signalHistos2[iecal][iclass]->Integral(1,signalHistos2[iecal][iclass]->GetNbinsX())/(float)signalHistos1[iecal][iclass]->Integral(1,signalHistos1[iecal][iclass]->GetNbinsX()));
    }

    double  max=TMath::Max(signalHistos1[iecal][0]->GetMaximum(),signalHistos1[iecal][1]->GetMaximum());
    max=TMath::Max(signalHistos2[iecal][0]->GetMaximum(),max);
    max=TMath::Max(signalHistos2[iecal][1]->GetMaximum(),max);
    signalHistos1[iecal][iclass]->SetMaximum(max+sqrt(max));

    double  min=TMath::Min(signalHistos1[iecal][0]->GetMinimum(),signalHistos1[iecal][1]->GetMinimum());
    signalHistos1[iecal][iclass]->SetMinimum(-2.);

    signalHistos1[iecal][iclass]->SetTitle("");
    signalHistos1[iecal][iclass]->GetXaxis()->SetTitle(axistitle);
    
    float nbins = signalHistos1[iecal][iclass]->GetNbinsX();
    float xmax = signalHistos1[iecal][iclass]->GetXaxis()->GetXmin();
    float xmin = signalHistos1[iecal][iclass]->GetXaxis()->GetXmax();
    char yaxDiv[200];
    sprintf(yaxDiv,"%.3f", fabs(xmax-xmin)/nbins);
     
    signalHistos1[iecal][iclass]->GetYaxis()->SetTitle("Events/"+TString(yaxDiv));
    }
    TString nameH0 = TString(signalHistos1[iecal][0]->GetName());
    TString nameH1 = TString(signalHistos1[iecal][1]->GetName());
    if(!nameH1.Contains("IPhiIPhi")) signalHistos1[iecal][1]->Draw("hist");
    if(!nameH0.Contains("fBrem")) {
      if(nameH1.Contains("IPhiIPhi")) signalHistos1[iecal][0]->Draw("hist");
      else signalHistos1[iecal][0]->Draw("hist same");
    }
    if(!nameH1.Contains("IPhiIPhi")) signalHistos2[iecal][1]->Draw("pE1 same");
    if(!nameH0.Contains("fBrem")) signalHistos2[iecal][0]->Draw("pE1 same");

    TPaveText pt1(0.1,0.92,0.8,0.97,"NDC");
    //	  pt1.SetTextFont(72);
    pt1.SetTextSize(0.05);
    pt1.SetTextAlign(12);
    pt1.SetFillColor(0);
    pt1.SetBorderSize(0);
    pt1.AddText("CMS Preliminary 2010, #sqrt{s}=7 TeV L_{int}="+intLumi+" pb^{-1}");
    //    pt1.Draw();
    //    c1.Update();

    TPaveText pt2(0.82,0.92,0.95,0.97,"NDC");
    pt2.SetTextSize(0.05);
    pt2.SetTextAlign(12);
    pt2.SetFillColor(0);
    pt2.SetBorderSize(0);
    if(iecal==0) pt2.AddText("Barrel");
    else pt2.AddText("Endcap");
    pt2.Draw();
    c1.Update();

    TLegend* leg = new TLegend(0.65,0.75,0.90,0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    if(MODE==0) {
      leg->AddEntry(signalHistos1[iecal][0],"MC 0 brem");
      leg->AddEntry(signalHistos1[iecal][1],"MC brem");
      leg->AddEntry(signalHistos2[iecal][0],"data 0 brem");
      leg->AddEntry(signalHistos2[iecal][1],"data brem");
    } else {
      leg->AddEntry(signalHistos1[iecal][0],"Signal MC 0 brem");
      leg->AddEntry(signalHistos1[iecal][1],"Signal MC brem");
      leg->AddEntry(signalHistos2[iecal][0],"Bkg MC 0 brem");
      leg->AddEntry(signalHistos2[iecal][1],"Bkg MC brem");
    }
    leg->Draw();

    char epsfilename[200];
    char rootfilename[200];
    char suffix[200];
    if(MODE==1) sprintf(suffix,"_sigVsBkg");
    else sprintf(suffix,"");
    if(iecal==0) {
      sprintf(epsfilename,"sPlotsData/%s_EB%s.png",namevar,suffix);
      sprintf(rootfilename,"sPlotsData/%s_EB%s.eps",namevar,suffix);
    }
    else if(iecal==1) {
      sprintf(epsfilename,"sPlotsData/%s_EE%s.png",namevar,suffix); 
      sprintf(rootfilename,"sPlotsData/%s_EE%s.eps",namevar,suffix);
    }
    c1.SaveAs(epsfilename);
    c1.SaveAs(rootfilename);

  }

}
