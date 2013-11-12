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
#include <TPaveText.h>
#include <TString.h>
#include <TAxis.h>

using namespace std;

void makePlots(TH1F *signalHistos1[2], TH1F *signalHistos2[2], const char *namevar, const char *axistitle);

int main(int argc, char* argv[]) {

  char inputFileNameSignal1[150];
  char inputFileNameSignal2[150];
  if ( argc < 3 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: MakeNotePdfPlots fileSignalMC.root fileSignalDATA.root" << std::endl;
    return 1;
  }
  strcpy(inputFileNameSignal1,argv[1]);
  strcpy(inputFileNameSignal2,argv[2]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *signalFile1 = TFile::Open(inputFileNameSignal1);
  signalFile1->cd();

  TH1F *dPhiClassEle1[2];
  TH1F *dEtaClassEle1[2];
  TH1F *EoPClassEle1[2];
  TH1F *HoEClassEle1[2];
  TH1F *sigmaIEtaIEtaClassEle1[2];
  TH1F *fbremClassEle1[2];
  TH1F *etaClassEle1[2];
  TH1F *phiClassEle1[2];
  TH1F *chargeClassEle1[2];
  TH1F *lhClassEle1[2];

  char histo[200];

  for(int iecal=0; iecal<2; iecal++) {

    sprintf(histo,"dPhiClass_electrons_%d",iecal);
    dPhiClassEle1[iecal]     = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"dEtaClass_electrons_%d",iecal);
    dEtaClassEle1[iecal]        = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"EoPClass_electrons_%d",iecal);
    EoPClassEle1[iecal]      = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"HoEClass_electrons_%d",iecal);
    HoEClassEle1[iecal]         = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
    sigmaIEtaIEtaClassEle1[iecal] = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"fbremClass_electrons_%d",iecal);
    fbremClassEle1[iecal] = (TH1F*)gDirectory->Get(histo); 
    // eta is the same for 0 and 1
    sprintf(histo,"etaClass_electrons");
    etaClassEle1[iecal]     = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"phiClass_electrons_%d",iecal);
    phiClassEle1[iecal] = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"chargeClass_electrons_%d",iecal);
    chargeClassEle1[iecal] = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"lhClass_electrons_%d",iecal);
    lhClassEle1[iecal] = (TH1F*)gDirectory->Get(histo); 
    
    dPhiClassEle1[iecal]   -> SetMinimum(0.001);
    dEtaClassEle1[iecal]   -> SetMinimum(0.001);
    EoPClassEle1[iecal]    -> SetMinimum(0.001);
    HoEClassEle1[iecal]   -> SetMinimum(0.001);
    sigmaIEtaIEtaClassEle1[iecal] -> SetMinimum(0.1);
    fbremClassEle1[iecal] -> SetMinimum(0.001);
    etaClassEle1[iecal] -> SetMinimum(0.001);
    phiClassEle1[iecal] -> SetMinimum(0.001);
    chargeClassEle1[iecal] -> SetMinimum(0.001);
    lhClassEle1[iecal] -> SetMinimum(0.001);
    
  }


  TFile *signalFile2 = TFile::Open(inputFileNameSignal2);
  signalFile2->cd();

  TH1F *dPhiClassEle2[2];
  TH1F *dEtaClassEle2[2];
  TH1F *EoPClassEle2[2];
  TH1F *HoEClassEle2[2];
  TH1F *sigmaIEtaIEtaClassEle2[2];
  TH1F *fbremClassEle2[2];
  TH1F *etaClassEle2[2];
  TH1F *phiClassEle2[2];
  TH1F *chargeClassEle2[2];
  TH1F *lhClassEle2[2];

  for(int iecal=0; iecal<2; iecal++) {
      
    sprintf(histo,"dPhiClass_electrons_%d",iecal);
    dPhiClassEle2[iecal]     = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"dEtaClass_electrons_%d",iecal);
    dEtaClassEle2[iecal]        = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"EoPClass_electrons_%d",iecal);
    EoPClassEle2[iecal]      = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"HoEClass_electrons_%d",iecal);
    HoEClassEle2[iecal]         = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
    sigmaIEtaIEtaClassEle2[iecal] = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"fbremClass_electrons_%d",iecal);
    fbremClassEle2[iecal] = (TH1F*)gDirectory->Get(histo); 
    // eta is the same for 0 and 1
    sprintf(histo,"etaClass_electrons",iecal);
    etaClassEle2[iecal]     = (TH1F*)gDirectory->Get(histo);
    sprintf(histo,"phiClass_electrons_%d",iecal);
    phiClassEle2[iecal] = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"chargeClass_electrons_%d",iecal);
    chargeClassEle2[iecal] = (TH1F*)gDirectory->Get(histo); 
    sprintf(histo,"lhClass_electrons_%d",iecal);
    lhClassEle2[iecal] = (TH1F*)gDirectory->Get(histo); 

    dPhiClassEle2[iecal]   -> SetMinimum(0.001);
    dEtaClassEle2[iecal]   -> SetMinimum(0.001);
    EoPClassEle2[iecal]    -> SetMinimum(0.001);
    HoEClassEle2[iecal]   -> SetMinimum(0.001);
    sigmaIEtaIEtaClassEle2[iecal] -> SetMinimum(0.1);
    fbremClassEle2[iecal]   -> SetMinimum(0.001);
    etaClassEle2[iecal] -> SetMinimum(0.001);
    phiClassEle2[iecal] -> SetMinimum(0.001);
    chargeClassEle2[iecal] -> SetMinimum(0.001);
    lhClassEle2[iecal] -> SetMinimum(0.001);

  }


  makePlots(dPhiClassEle1,dPhiClassEle2,"dPhiIn","#Delta #phi_{in} (rad)");
  makePlots(dEtaClassEle1,dEtaClassEle2,"dEtaIn","#Delta #eta_{in} (rad)");
  makePlots(EoPClassEle1,EoPClassEle2,"EoP","E_{SC}/p_{in}");
  makePlots(HoEClassEle1,HoEClassEle2,"HoE","H/E");
  makePlots(sigmaIEtaIEtaClassEle1,sigmaIEtaIEtaClassEle2,"sigmaIEtaIEta","#sigma_{#eta#eta} (rad)");
  makePlots(fbremClassEle1,fbremClassEle2,"fbrem","f_{brem}");
  makePlots(etaClassEle1,etaClassEle2,"eta","#eta");
  //  makePlots(phiClassEle1,phiClassEle2,"phi","#phi");
  makePlots(chargeClassEle1,chargeClassEle2,"charge","charge");
  makePlots(lhClassEle1,lhClassEle2,"lh","LH output");

}

void makePlots(TH1F *signalHistos1[2], TH1F *signalHistos2[2], const char *namevar, const char *axistitle) {

  TString intLumi="2.88";

  int iptbin=1;

  for(int iecal=0; iecal<2; iecal++) {

    TCanvas c1;
    // c1.SetLogy();

    // MC
    if(signalHistos1[iecal]) {
      signalHistos1[iecal]->SetLineWidth(2);
      signalHistos1[iecal]->SetLineColor(kRed+1);
      signalHistos1[iecal]->SetFillColor(kYellow+1);
//       signalHistos1[iecal]->SetLineColor(kViolet+3);
//       signalHistos1[iecal]->SetFillColor(kViolet);
//      signalHistos1[iecal]->SetLineColor(kAzure+3);
//      signalHistos1[iecal]->SetFillColor(kAzure+6);
      signalHistos1[iecal]->Scale((float)signalHistos2[iecal]->Integral(1,signalHistos2[iecal]->GetNbinsX())/(float)signalHistos1[iecal]->Integral(1,signalHistos1[iecal]->GetNbinsX()));
      std::cout << "Normalized to " << (float)signalHistos2[iecal]->Integral(1,signalHistos2[iecal]->GetNbinsX()) << std::endl;
    }
    // data
    if(signalHistos2[iecal]) {
      signalHistos2[iecal]->SetLineWidth(2);
      signalHistos2[iecal]->SetLineColor(kBlack);
      signalHistos2[iecal]->SetMarkerColor(kBlack);
      signalHistos2[iecal]->SetMarkerStyle(8);
      signalHistos2[iecal]->SetMarkerSize(1.5);
      // signalHistos2[iecal]->Scale(1.0/signalHistos2[iecal]->Integral());
    }

    double  max=TMath::Max(signalHistos1[iecal]->GetMaximum(),signalHistos2[iecal]->GetMaximum());
    signalHistos1[iecal]->SetMaximum(max+sqrt(max));

    double  min=TMath::Min(signalHistos1[iecal]->GetMinimum(),signalHistos2[iecal]->GetMinimum());
    signalHistos1[iecal]->SetMinimum(-2.);

    signalHistos1[iecal]->SetTitle("");
    signalHistos1[iecal]->GetXaxis()->SetTitle(axistitle);
    
    float nbins = signalHistos1[iecal]->GetNbinsX();
    float xmax = signalHistos1[iecal]->GetXaxis()->GetXmin();
    float xmin = signalHistos1[iecal]->GetXaxis()->GetXmax();
    char yaxDiv[200];
    sprintf(yaxDiv,"%.3f", fabs(xmax-xmin)/nbins);
     
    signalHistos1[iecal]->GetYaxis()->SetTitle("Events/"+TString(yaxDiv));
    signalHistos1[iecal]->Draw("hist");
    signalHistos2[iecal]->Draw("pE1 same");

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
    //    leg->AddEntry(signalHistos1[iecal],"MC");
    //    leg->AddEntry(signalHistos2[iecal],"data");
    leg->AddEntry(signalHistos1[iecal],"single e, p_{T}>20 GeV");
    leg->AddEntry(signalHistos2[iecal],"W #rightarrow e#nu, p_{T}>20 GeV");
    leg->Draw();

    char epsfilename[200];
    char rootfilename[200];
    if(iecal==0) {
      sprintf(epsfilename,"sPlotsData/%s_EB.png",namevar);
      sprintf(rootfilename,"sPlotsData/%s_EB.eps",namevar);
    }
    else if(iecal==1) {
      sprintf(epsfilename,"sPlotsData/%s_EE.png",namevar); 
      sprintf(rootfilename,"sPlotsData/%s_EE.eps",namevar);
    }
    c1.SaveAs(epsfilename);
    c1.SaveAs(rootfilename);

  }

}
