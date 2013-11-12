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
#include <TLegend.h>
#include <TCanvas.h>

// Offline analysis includes
#include "CommonTools/include/EfficiencyEvaluator.hh"

#include <string>
#include <sstream>
#include <fstream>

using namespace std;

std::vector<TString> EgammaCutBasedIDWPs;
//std::vector<TString> EgammaCiCBasedIDWPs;
std::vector<TString> EgammaLHBasedIDWPs;

void drawEta(const char* filename, TString type);
void drawPt(const char* filename, TString type);
void drawPU(const char* filename, TString type);

int main(int argc, char* argv[]) {

  TString IDpart;
   if (argv[4])
     {
       char idpar[100];
       strcpy(idpar,argv[4]);
       IDpart=TString(idpar);
     }
   else
     IDpart="";

   EgammaCutBasedIDWPs.push_back("WP95"+TString(IDpart));;
   EgammaCutBasedIDWPs.push_back("WP90"+TString(IDpart));;
   EgammaCutBasedIDWPs.push_back("WP85"+TString(IDpart));;
  EgammaCutBasedIDWPs.push_back("WP80"+TString(IDpart));;
  EgammaCutBasedIDWPs.push_back("WP70"+TString(IDpart));;    
  //  EgammaCiCBasedIDWPs.push_back(TString("CiCVeryLoose")+TString(IDpart));;
  //  EgammaCiCBasedIDWPs.push_back(TString("CiCLoose")+TString(IDpart));;
  //EgammaCiCBasedIDWPs.push_back(TString("CiCLoose")+TString(IDpart));;
  // EgammaCiCBasedIDWPs.push_back(TString("CiCTight")+TString(IDpart));;
  //EgammaCiCBasedIDWPs.push_back(TString("CiCSuperTight")+TString(IDpart));;
  //EgammaCiCBasedIDWPs.push_back(TString("CiCHyperTight")+TString(IDpart));;
  //  EgammaCiCBasedIDWPs.push_back(TString("CiCHyperTight2")+TString(IDpart));;

  EgammaLHBasedIDWPs.push_back(TString("LHVeryLoose")+TString(IDpart));;
  EgammaLHBasedIDWPs.push_back(TString("LHLoose")+TString(IDpart));;
  EgammaLHBasedIDWPs.push_back(TString("LHMedium")+TString(IDpart));;
  EgammaLHBasedIDWPs.push_back(TString("LHTight")+TString(IDpart));;
  EgammaLHBasedIDWPs.push_back(TString("LHHyperTight")+TString(IDpart));;

  char inputFileNameEta[150];
  char inputFileNamePt[150];
  char inputFileNamePU[150];
  char type[150];
  if ( argc < 4 ){
    std::cout << "missing argument!" << std::endl; 
    std::cout << "usage: compareEff inputFileEtaPrefix inputFilePtPrefix inputFilePUPrefix <optional id part>" << std::endl;
    return 1;
  }
  strcpy(inputFileNameEta,argv[1]);
  strcpy(inputFileNamePt,argv[2]);
  strcpy(inputFileNamePU,argv[3]);

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  char input[300];
//   sprintf(input,"%s.root",inputFileNameEta);
//   drawEta(input,"");
  sprintf(input,"%sHighPt.root",inputFileNameEta);
  drawEta(input,"HighPt");
  sprintf(input,"%sLowPt.root",inputFileNameEta);
  drawEta(input,"LowPt");
//   sprintf(input,"%s.root",inputFileNamePt);
//   drawPt(input, "");
  sprintf(input,"%sBarrel.root",inputFileNamePt);
  drawPt(input, "Barrel");
  sprintf(input,"%sEndcap.root",inputFileNamePt);
  drawPt(input, "Endcap");

//   sprintf(input,"%s.root",inputFileNamePU);
//   drawPU(input,"");
  sprintf(input,"%sHighPt.root",inputFileNamePU);
  drawPU(input,"HighPt");
  sprintf(input,"%sLowPt.root",inputFileNamePU);
  drawPU(input,"LowPt");

}


void drawEta(const char* filename,TString type) {

  TFile *efficiencyFileEta = TFile::Open(filename);

  //  TH1F *GenEta = (TH1F*)efficiencyFileEta->Get("GenEta"+type);
  TH1F *RecoEta = (TH1F*)efficiencyFileEta->Get("RecoEta"+type);

  std::vector<TH1F*> CutIdEta;
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"Eta"+type);
      CutIdEta.push_back(aHisto);
//       aHisto = (TH1F*)efficiencyFileEta->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDEta"+type);
//       CutIdEta.push_back(aHisto);
    }
  std::vector<TH1F*> LHIdEta;
  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"Eta"+type);
      LHIdEta.push_back(aHisto);
//       aHisto = (TH1F*)efficiencyFileEta->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDEta"+type);
//       LHIdEta.push_back(aHisto);
    }
//   std::vector<TH1F*> CiCIdEta;
//   for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
//     {
//       TH1F* aHisto = (TH1F*)efficiencyFileEta->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Eta"+type);
//       CiCIdEta.push_back(aHisto);
// //       aHisto = (TH1F*)efficiencyFileEta->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDEta"+type);
// //       CiCIdEta.push_back(aHisto);
//     }


  EfficiencyEvaluator ElectronEfficiencyEta("eleeff-eta_"+type+".root");
  //  ElectronEfficiencyEta.AddNumerator(GenEta);
  //  ElectronEfficiencyEta.AddNumerator(RecoEta);
  for (int i=0;i<CutIdEta.size();++i)
    ElectronEfficiencyEta.AddNumerator(CutIdEta[i]);
  for (int i=0;i<LHIdEta.size();++i)
    ElectronEfficiencyEta.AddNumerator(LHIdEta[i]);
//   for (int i=0;i<CiCIdEta.size();++i)
//     ElectronEfficiencyEta.AddNumerator(CiCIdEta[i]);
  ElectronEfficiencyEta.SetDenominator(RecoEta);
  ElectronEfficiencyEta.ComputeEfficiencies();
  
  vector<TH1F*> efficiency = ElectronEfficiencyEta.GetCumulativeEfficiencies();
  vector<float> effAverageEB = ElectronEfficiencyEta.GetCumulativeEfficiencyAverages(13,47);
  vector<float> effAverageEEm = ElectronEfficiencyEta.GetCumulativeEfficiencyAverages(1,12);
  vector<float> effAverageEEp = ElectronEfficiencyEta.GetCumulativeEfficiencyAverages(48,60);

  std::ostringstream effCutEB;
  std::ostringstream effCutEE;
  std::ostringstream effLHEB;
  std::ostringstream effLHEE;

  effCutEB << "float eff_EB_CutId_" << type << "[" << CutIdEta.size() << "] = { ";
  effCutEE << "float eff_EE_CutId_" << type << "[" << CutIdEta.size() << "] = { ";

  effLHEB << "float eff_EB_LHId_" << type << "[" << CutIdEta.size() << "] = { ";
  effLHEE << "float eff_EE_LHId_" << type << "[" << CutIdEta.size() << "] = { ";

  for (int tightness=0;tightness<CutIdEta.size();++tightness)
    { 
      effCutEB << effAverageEB[0+tightness] ;
      effCutEE << (effAverageEEm[0+tightness]+effAverageEEp[0+tightness])/2 ;
      effLHEB << effAverageEB[CutIdEta.size()+tightness] ;
      effLHEE << (effAverageEEm[CutIdEta.size()+tightness]+effAverageEEp[CutIdEta.size()+tightness])/2. ;
      if (tightness==CutIdEta.size()-1)
	{
	  effCutEB << " };"  ;
	  effCutEE << " };"  ;
	  effLHEB << " };"  ;
	  effLHEE << " };"  ;
	}
      else
	{
	  effCutEB << " , " ;
	  effCutEE << " , " ;
	  effLHEB << " , " ;
	  effLHEE << " , " ;
	}
    }

  std::ostringstream fileName;
  fileName << "efficiencies-" << type << ".txt";

  std::ofstream out;
  out.open(fileName.str().c_str());
  out << "//************* EB " << type << " efficiencies *************" << std::endl;
  out << effCutEB.str() << std::endl;
  out << effLHEB.str() << std::endl;

  out << "//************* EE " << type << " efficiencies *************" << std::endl;
  out << effCutEE.str() << std::endl;
  out << effLHEE.str() << std::endl;
  

  for (int tightness=0;tightness<CutIdEta.size();tightness++)
    { 
      TCanvas c1;
      efficiency[0+tightness]->SetLineColor(1);
      efficiency[0+tightness]->SetMarkerColor(1);
      efficiency[CutIdEta.size()+tightness]->SetLineColor(2);
      efficiency[CutIdEta.size()+tightness]->SetMarkerColor(2);
//       efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->SetLineColor(4);
//       efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->SetMarkerColor(4);
      //  efficiency[6]->SetLineColor(8);

      char effCutIdEB[50];
      sprintf(effCutIdEB," (#epsilon EB = %.2f) ", effAverageEB[0+tightness]);
      char effCutIdEE[50];
      sprintf(effCutIdEE," (#epsilon EE = %.2f) ", (effAverageEEm[0+tightness]+effAverageEEp[0+tightness])/2.);
      char effLHIdEB[50];
      sprintf(effLHIdEB," (#epsilon EB = %.2f) ", effAverageEB[CutIdEta.size()+tightness]);
      char effLHIdEE[50];
      sprintf(effLHIdEE," (#epsilon EE = %.2f) ", (effAverageEEm[CutIdEta.size()+tightness]+effAverageEEp[CutIdEta.size()+tightness])/2.);
      
      TLegend* leg = new TLegend(0.35,0.25,0.60,0.50);
      leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
      leg->SetFillColor(0);
      leg->AddEntry(efficiency[0+tightness],(EgammaCutBasedIDWPs[0+tightness]+TString(effCutIdEB)).Data());
      leg->AddEntry(efficiency[0+tightness],(EgammaCutBasedIDWPs[0+tightness]+TString(effCutIdEE)).Data());
      leg->AddEntry(efficiency[CutIdEta.size()+tightness],(EgammaLHBasedIDWPs[0+tightness]+TString(effLHIdEB)).Data());
      leg->AddEntry(efficiency[CutIdEta.size()+tightness],(EgammaLHBasedIDWPs[0+tightness]+TString(effLHIdEE)).Data());
      //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
      //  leg->AddEntry(efficiency[CutIdEta.size()+LHIdEta.size()+tightness],EgammaCiCBasedIDWPs[0+tightness].Data());
      
      efficiency[0+tightness]->SetMinimum(0.1);
      efficiency[0+tightness]->SetMaximum(1.0);
      efficiency[0+tightness]->SetMarkerStyle(20);
      efficiency[0+tightness]->SetMarkerSize(1.05);
      efficiency[0+tightness]->SetTitle("");
      efficiency[0+tightness]->GetXaxis()->SetTitle("electron #eta");
      efficiency[0+tightness]->GetYaxis()->SetTitle("Efficiency");
      efficiency[0+tightness]->Draw("PEhist");
      efficiency[CutIdEta.size()+tightness]->Draw("PEhistsame");
      efficiency[CutIdEta.size()+tightness]->SetMarkerStyle(21);
      efficiency[CutIdEta.size()+tightness]->SetMarkerSize(1.05);
      //  efficiency[3+tightness]->Draw("same hist");
//       efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->Draw("PEhistsame");
//       efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->SetMarkerStyle(21);
//       efficiency[CutIdEta.size()+LHIdEta.size()+tightness]->SetMarkerSize(1.05);
      leg->Draw();
      
      char tightLevel[10];
      sprintf(tightLevel,"%d",tightness);
      c1.SaveAs("eleeff-eta-"+type+"-tightness"+TString(tightLevel)+".eps");
      c1.SaveAs("eleeff-eta-"+type+"-tightness"+TString(tightLevel)+".root");
      c1.SaveAs("eleeff-eta-"+type+"-tightness"+TString(tightLevel)+".png");
    }

}

void drawPt(const char* filename, TString type) {
  TFile *efficiencyFilePt = TFile::Open(filename);

  //  TH1F *GenPt = (TH1F*)efficiencyFilePt->Get("GenPt"+type);
  TH1F *RecoPt = (TH1F*)efficiencyFilePt->Get("RecoPt"+type);

  std::vector<TH1F*> CutIdPt;
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"Pt"+type);
      CutIdPt.push_back(aHisto);
    }
  std::vector<TH1F*> LHIdPt;
  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"Pt"+type);
      LHIdPt.push_back(aHisto);
    }
//   std::vector<TH1F*> CiCIdPt;
//   for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
//     {
//       TH1F* aHisto = (TH1F*)efficiencyFilePt->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"Pt"+type);
//       CiCIdPt.push_back(aHisto);
//     }




  EfficiencyEvaluator ElectronEfficiencyPt("eleeff-pt"+type+".root");
  //  ElectronEfficiencyPt.AddNumerator(GenPt);
  //  ElectronEfficiencyPt.AddNumerator(RecoPt);
  for (int i=0;i<CutIdPt.size();++i)
    ElectronEfficiencyPt.AddNumerator(CutIdPt[i]);
  for (int i=0;i<LHIdPt.size();++i)
    ElectronEfficiencyPt.AddNumerator(LHIdPt[i]);
//   for (int i=0;i<CiCIdPt.size();++i)
//     ElectronEfficiencyPt.AddNumerator(CiCIdPt[i]);
  ElectronEfficiencyPt.SetDenominator(RecoPt);
  ElectronEfficiencyPt.ComputeEfficiencies();
  
  vector<TH1F*> efficiency = ElectronEfficiencyPt.GetCumulativeEfficiencies();
  vector<float> effAverage = ElectronEfficiencyPt.GetCumulativeEfficiencyAverages();
      
  for (int tightness=0;tightness<CutIdPt.size();++tightness)
    { 
      TCanvas c1;
      efficiency[0+tightness]->SetLineColor(1);
      efficiency[0+tightness]->SetMarkerColor(1);
      efficiency[CutIdPt.size()+tightness]->SetLineColor(2);
      efficiency[CutIdPt.size()+tightness]->SetMarkerColor(2);
//       efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->SetLineColor(4);
//       efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->SetMarkerColor(4);
      //  efficiency[6]->SetLineColor(8);
      
      char effCutId[50];
      sprintf(effCutId," (#epsilon = %.2f) ", effAverage[0+tightness]);
      char effLHId[50];
      sprintf(effLHId," (#epsilon = %.2f) ", effAverage[CutIdPt.size()+tightness]);

      TLegend* leg = new TLegend(0.35,0.25,0.60,0.50);
      leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
      leg->SetFillColor(0);
      leg->AddEntry(efficiency[0+tightness],(EgammaCutBasedIDWPs[0+tightness]+TString(effCutId)).Data());
      leg->AddEntry(efficiency[CutIdPt.size()+tightness],(EgammaLHBasedIDWPs[0+tightness]+TString(effLHId)).Data());
      //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
      //  leg->AddEntry(efficiency[CutIdPt.size()+LHIdPt.size()+tightness],EgammaCiCBasedIDWPs[0+tightness].Data());
      
      efficiency[0+tightness]->SetMinimum(0.1);
      efficiency[0+tightness]->SetMaximum(1.0);
      efficiency[0+tightness]->SetMarkerStyle(20);
      efficiency[0+tightness]->SetMarkerSize(1.05);
      efficiency[0+tightness]->SetTitle("");
      efficiency[0+tightness]->GetXaxis()->SetTitle("electron #pt");
      efficiency[0+tightness]->GetYaxis()->SetTitle("Efficiency");
      efficiency[0+tightness]->Draw("PEhist");
      efficiency[CutIdPt.size()+tightness]->Draw("PEhistsame");
      efficiency[CutIdPt.size()+tightness]->SetMarkerStyle(21);
      efficiency[CutIdPt.size()+tightness]->SetMarkerSize(1.05);
      //  efficiency[3+tightness]->Draw("same hist");
//       efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->Draw("PEhistsame");
//       efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->SetMarkerStyle(21);
//       efficiency[CutIdPt.size()+LHIdPt.size()+tightness]->SetMarkerSize(1.05);
      leg->Draw();
      
      char tightLevel[10];
      sprintf(tightLevel,"%d",tightness);
      c1.SaveAs("eleeff-pt-"+type+"-tightness"+TString(tightLevel)+".eps");
      c1.SaveAs("eleeff-pt-"+type+"-tightness"+TString(tightLevel)+".root");
      c1.SaveAs("eleeff-pt-"+type+"-tightness"+TString(tightLevel)+".png");
    }


}

void drawPU(const char* filename,TString type) {

  TFile *efficiencyFilePU = TFile::Open(filename);

  //  TH1F *GenPU = (TH1F*)efficiencyFilePU->Get("GenPU"+type);
  TH1F *RecoPU = (TH1F*)efficiencyFilePU->Get("RecoPU"+type);

  std::vector<TH1F*> CutIdPU;
  for (int i=0;i<EgammaCutBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePU->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"PU"+type);
      CutIdPU.push_back(aHisto);
//       aHisto = (TH1F*)efficiencyFilePU->Get("CutId"+TString(EgammaCutBasedIDWPs[i])+"OnlyIDPU"+type);
//       CutIdPU.push_back(aHisto);
    }
  std::vector<TH1F*> LHIdPU;
  for (int i=0;i<EgammaLHBasedIDWPs.size();++i)
    {
      TH1F* aHisto = (TH1F*)efficiencyFilePU->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"PU"+type);
      LHIdPU.push_back(aHisto);
//       aHisto = (TH1F*)efficiencyFilePU->Get("LHId"+TString(EgammaLHBasedIDWPs[i])+"OnlyIDPU"+type);
//       LHIdPU.push_back(aHisto);
    }
//   std::vector<TH1F*> CiCIdPU;
//   for (int i=0;i<EgammaCiCBasedIDWPs.size();++i)
//     {
//       TH1F* aHisto = (TH1F*)efficiencyFilePU->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"PU"+type);
//       CiCIdPU.push_back(aHisto);
// //       aHisto = (TH1F*)efficiencyFilePU->Get("CiCId"+TString(EgammaCiCBasedIDWPs[i])+"OnlyIDPU"+type);
// //       CiCIdPU.push_back(aHisto);
//     }




  EfficiencyEvaluator ElectronEfficiencyPU("eleeff-PU_"+type+".root");
  //  ElectronEfficiencyPU.AddNumerator(GenPU);
  //  ElectronEfficiencyPU.AddNumerator(RecoPU);
  for (int i=0;i<CutIdPU.size();++i)
    ElectronEfficiencyPU.AddNumerator(CutIdPU[i]);
  for (int i=0;i<LHIdPU.size();++i)
    ElectronEfficiencyPU.AddNumerator(LHIdPU[i]);
//   for (int i=0;i<CiCIdPU.size();++i)
//     ElectronEfficiencyPU.AddNumerator(CiCIdPU[i]);
  ElectronEfficiencyPU.SetDenominator(RecoPU);
  ElectronEfficiencyPU.ComputeEfficiencies();
  
  vector<TH1F*> efficiency = ElectronEfficiencyPU.GetCumulativeEfficiencies();
  vector<float> effAverage = ElectronEfficiencyPU.GetCumulativeEfficiencyAverages();

  for (int tightness=0;tightness<CutIdPU.size();tightness++)
    { 
      TCanvas c1;
      efficiency[0+tightness]->SetLineColor(1);
      efficiency[0+tightness]->SetMarkerColor(1);
      efficiency[CutIdPU.size()+tightness]->SetLineColor(2);
      efficiency[CutIdPU.size()+tightness]->SetMarkerColor(2);
//       efficiency[CutIdPU.size()+LHIdPU.size()+tightness]->SetLineColor(4);
//       efficiency[CutIdPU.size()+LHIdPU.size()+tightness]->SetMarkerColor(4);
      //  efficiency[6]->SetLineColor(8);
     
      char effCutId[50];
      sprintf(effCutId," (#epsilon = %.2f) ", effAverage[0+tightness]);
      char effLHId[50];
      sprintf(effLHId," (#epsilon = %.2f) ", effAverage[CutIdPU.size()+tightness]);
 
      TLegend* leg = new TLegend(0.35,0.25,0.60,0.50);
      leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
      leg->SetFillColor(0);
      leg->AddEntry(efficiency[0+tightness],(EgammaCutBasedIDWPs[0+tightness]+TString(effCutId)).Data());
      leg->AddEntry(efficiency[CutIdPU.size()+tightness],(EgammaLHBasedIDWPs[0+tightness]+TString(effLHId)).Data());
      //  leg->AddEntry(efficiency[3],"Tight CIC Cuts");
      //  leg->AddEntry(efficiency[CutIdPU.size()+LHIdPU.size()+tightness],EgammaCiCBasedIDWPs[0+tightness].Data());
      
      efficiency[0+tightness]->SetMinimum(0.1);
      efficiency[0+tightness]->SetMaximum(1.0);
      efficiency[0+tightness]->SetMarkerStyle(20);
      efficiency[0+tightness]->SetMarkerSize(1.05);
      efficiency[0+tightness]->SetTitle("");
      efficiency[0+tightness]->GetXaxis()->SetTitle("n PU");
      efficiency[0+tightness]->GetYaxis()->SetTitle("Efficiency");
      efficiency[0+tightness]->Draw("PEhist");
      efficiency[CutIdPU.size()+tightness]->Draw("PEhistsame");
      efficiency[CutIdPU.size()+tightness]->SetMarkerStyle(21);
      efficiency[CutIdPU.size()+tightness]->SetMarkerSize(1.05);
      //  efficiency[3+tightness]->Draw("same hist");
//       efficiency[CutIdPU.size()+LHIdPU.size()+tightness]->Draw("PEhistsame");
//       efficiency[CutIdPU.size()+LHIdPU.size()+tightness]->SetMarkerStyle(21);
//       efficiency[CutIdPU.size()+LHIdPU.size()+tightness]->SetMarkerSize(1.05);
      leg->Draw();
      
      char tightLevel[10];
      sprintf(tightLevel,"%d",tightness);
      c1.SaveAs("eleeff-PU-"+type+"-tightness"+TString(tightLevel)+".eps");
      c1.SaveAs("eleeff-PU-"+type+"-tightness"+TString(tightLevel)+".root");
      c1.SaveAs("eleeff-PU-"+type+"-tightness"+TString(tightLevel)+".png");
    }

}
