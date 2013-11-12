#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TH2F.h"
#include "TString.h"
#include <vector>
#include <iostream>

void makeEffScaleFact() {

  gStyle->SetPalette(1);

  // data contains TH2F with fit results
  TFile *fileData = TFile::Open("egammaresults/tpfits/hzzidref/2011/newhzzWP_eff_Data7TeV_Kine_Pt10To1000.root");
  TCanvas *cData = (TCanvas*)fileData->Get("eleIDdir/etapt/fit_eff_plots/pt_abseta_PLOT_vertices_bin0");
  TH2F *hData = (TH2F*)cData->GetPrimitive("pt_abseta_PLOT_vertices_bin0");
  std::cout << "Data histogram = " << hData << std::endl;

  TFile *fileMC = TFile::Open("egammaresults/tpfits/hzzidref/2011/newhzzWP_eff_MC7TeV_Kine_Pt10To1000.root");
  TCanvas *cMC = (TCanvas*)fileMC->Get("eleIDdir/etapt/cnt_eff_plots/pt_abseta_PLOT_vertices_bin0");
  TH2F *hMC = (TH2F*)cMC->GetPrimitive("pt_abseta_PLOT_vertices_bin0");
  std::cout << "MC histogram = " << hMC << std::endl;

  hData->Divide(hMC);

  TCanvas scaleF("scaleF","scaleF",1600,1067);
  hData->Draw();
  //  hData->Draw("text same");
  
  TFile *fileOut = TFile::Open("effsf_hzzref_2011data_10To1000.root","recreate");
  hData->Write();
  fileOut->Close();

}
