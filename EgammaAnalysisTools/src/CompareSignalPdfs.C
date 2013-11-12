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
#include <TTree.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TPad.h>
#include <RooDataSet.h>

char inputFileNameData[150];
char inputFileNameSigMC[150];
char inputFileNameQCDMC[150];

using namespace std;

void makePlots(const char *namevar, const char *axistitle, int nbins, float min, float max);

int ComparePdfs() {

  sprintf(inputFileNameData,"results/datasets/sPlots/Wenu.root");
  sprintf(inputFileNameSigMC,"results/isolation_trees/WJetsMADGRAPH_IsolVtx.root");
  sprintf(inputFileNameQCDMC,"results/isolation_trees/QCD_IsolVtx.root");

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  //[iecal:0=EB;1=EE][isample:0=sig,1=qcd][imc:0=mc,1=data]
  int nbins = 100;
  float deltaPhiMin  = -0.1;
  float deltaPhiMax  =  0.1;
  float deltaEtaMin     = -0.02;
  float deltaEtaMax     =  0.02;
  float EoPMin      =  0.0;
  float EoPMax      =  8.0;
  float HoEMin      =  0.0;
  float HoEMax      =  0.2;
  float sigmaIEtaIEtaMin = 0.0;
  float sigmaIEtaIEtaMax = 0.05;

  makePlots("deltaPhi","#Delta #phi_{in} (rad)",nbins,deltaPhiMin,deltaPhiMax);
  makePlots("deltaEta","#Delta #eta_{in} (rad)",nbins,deltaEtaMin,deltaEtaMax);
  makePlots("eop","E_{SC}/P_{in}",nbins,EoPMin,EoPMax);
  makePlots("see","#sigma_{i#eta i#eta}",nbins,sigmaIEtaIEtaMin,sigmaIEtaIEtaMax);

}

void makePlots(const char *namevar, const char *axistitle, int nbins, float min, float max) {

  TFile *dataFile = TFile::Open(inputFileNameData);
  RooDataSet* dataset = (RooDataSet*)dataFile->Get("dataset");

  TFile *sigFile = TFile::Open(inputFileNameSigMC);
  TTree *sigTree = (TTree*)sigFile->Get("T1");

  TFile *qcdFile = TFile::Open(inputFileNameQCDMC);
  TTree *qcdTree = (TTree*)qcdFile->Get("T1");

  TTree *mcTree = 0;

  TCanvas c1;
  // c1.SetLogy();

  TPad *npad = 0;
  
  TH1F *pdf[2][2][2];
  for(int iecal=0; iecal<2; iecal++) {
    for(int isample=0; isample<2; isample++) {
      for(int imc=0; imc<2; imc++) {
        
        char histo[200];
        sprintf(histo,"%s_%d_%d_%d",namevar,iecal,isample,imc);
        pdf[iecal][isample][imc]    = new TH1F(histo, histo, nbins, min, max);
        pdf[iecal][isample][imc]->Sumw2();
      }
    }
  }  

  for(int iecal=0; iecal<2; iecal++) {
          
    std::string cut;
    if(iecal==0) cut="(abs(eta)<1.479)";
    else cut="(abs(eta)>1.479)";

    for(int isample=0; isample<2; isample++) {

      std::string dataWeight;
      if(isample==0) {
        dataWeight = cut+std::string("*0.01*N_sig_sw");
        mcTree = sigTree;
      } else {
        dataWeight = cut+std::string("*0.01*N_qcd_sw");
        mcTree = qcdTree;
      }

      // fill MC
      TString name = pdf[iecal][isample][0]->GetName();
      mcTree->Project(name.Data(),namevar,cut.c_str());

      // fill data s-weighted
      name = pdf[iecal][isample][1]->GetName();
      TTree *treeData = (TTree*)dataset->tree();
      treeData->Project(name.Data(),namevar,dataWeight.c_str());
    
      // style
      pdf[iecal][isample][0]->SetLineWidth(2);
      pdf[iecal][isample][0]->SetLineColor(kRed+1);
      pdf[iecal][isample][0]->Scale(pdf[iecal][isample][1]->Integral()/pdf[iecal][isample][0]->Integral());
     
      float max = TMath::Max(pdf[iecal][isample][0]->GetMaximum(),pdf[iecal][isample][1]->GetMaximum());
      pdf[iecal][isample][0]->SetMaximum(max+sqrt(max));

      pdf[iecal][isample][1]->SetLineWidth(2);
      pdf[iecal][isample][1]->SetMarkerStyle(4);
      pdf[iecal][isample][1]->SetMarkerColor(kBlack);
      
      pdf[iecal][isample][0]->SetTitle("");
      pdf[iecal][isample][0]->GetXaxis()->SetTitle(axistitle);
      pdf[iecal][isample][0]->GetYaxis()->SetTitle("normalized to 10/nb");

    }

    c1.cd();
    // large plot for the signal
    pdf[iecal][0][0]->Draw("hist");
    pdf[iecal][0][1]->Draw("same pe1");

    // inset for the background
    npad = new TPad("npad", "", 0.6, 0.6, 0.9, 0.88);
    npad->Draw();
    npad->cd();
    // -- make enough space so that axis titles are not cut off
    npad->SetBottomMargin(0.25);
    npad->SetLeftMargin(0.25);

    pdf[iecal][1][0]->Draw("hist");
    pdf[iecal][1][1]->Draw("same pe1");

    // -- go back to the big canvas and put a legend there
    c1.cd();

    TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->AddEntry(pdf[iecal][0][0],"MC");
    leg->AddEntry(pdf[iecal][0][1],"data");
    leg->Draw();

    char epsfilename[200];
    char rootfilename[200];
    if(iecal==0) {
      sprintf(epsfilename,"%s_EB.eps",namevar);
      sprintf(rootfilename,"%s_EB.root",namevar);
    } else {
      sprintf(epsfilename,"%s_EE.eps",namevar); 
      sprintf(rootfilename,"%s_EE.root",namevar);
    }
    c1.SaveAs(epsfilename);
    c1.SaveAs(rootfilename);

  }

}
