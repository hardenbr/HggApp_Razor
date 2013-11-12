#include "TH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>

void efficiency(std::vector<TH1F*> histos_MC, std::vector<TH1F*> histos_data, int iecal);

int EffTest(int iecal) {

  TFile *fileData = TFile::Open("pdfs_histograms_data_x2.root");
  TFile *fileMC = TFile::Open("pdfs_histograms_MC_x2.root");
  
  std::vector<TH1F*> histos_MC;
  std::vector<TH1F*> histos_data;
  histos_MC.clear();
  histos_data.clear();

  char histo[200];
  sprintf(histo,"dPhiClass_electrons_%d",iecal);
  histos_MC.push_back((TH1F*)fileMC->Get(histo));
  histos_data.push_back((TH1F*)fileData->Get(histo));
  
  sprintf(histo,"dEtaClass_electrons_%d",iecal);
  histos_MC.push_back((TH1F*)fileMC->Get(histo));
  histos_data.push_back((TH1F*)fileData->Get(histo));
  
  sprintf(histo,"HoEClass_electrons_%d",iecal);
  histos_MC.push_back((TH1F*)fileMC->Get(histo));
  histos_data.push_back((TH1F*)fileData->Get(histo));
  
  sprintf(histo,"EoPClass_electrons_%d",iecal);
  histos_MC.push_back((TH1F*)fileMC->Get(histo));
  histos_data.push_back((TH1F*)fileData->Get(histo));
  
  sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
  histos_MC.push_back((TH1F*)fileMC->Get(histo));
  histos_data.push_back((TH1F*)fileData->Get(histo));
  
  efficiency(histos_MC, histos_data, iecal);    

  return 0;

}

void efficiency(std::vector<TH1F*> histos_MC, std::vector<TH1F*> histos_data, int iecal) {


  float start, stop, step;
  bool abs;
  int nsteps;

  TGraphErrors g_eff_MC, g_eff_Data, g_eff_residual;

  for(int i=0; i<(int)histos_MC.size(); ++i) {
    
    TString name(histos_data[i]->GetName());
    if(name.Contains("dPhi")) {
      start = 0.0;
      stop = 0.04;
      nsteps = 10;
      step = (stop-start)/float(nsteps);
      abs = true;
    } else if(name.Contains("dEta")) {
      start = 0.0;
      stop = 0.004;
      nsteps = 10;
      step = (stop-start)/float(nsteps);
      abs = true;
    } else if(name.Contains("sigma")) {
      if(iecal==0) {
        start = 0.005;
        stop = 0.012;
      } else {
        start = 0.02;
        stop = 0.03;
      }
      nsteps = 10;
      step = (stop-start)/float(nsteps);
      abs = false;
    } else if(name.Contains("HoE")) {
      if(iecal==0) {
        start = 0.0;
        stop = 0.07;
      } else {
        start = 0.0;
        stop = 0.06;
      }
      nsteps = 10;
      step = (stop-start)/float(nsteps);
      abs = false;
    } else if(name.Contains("EoP")) {
      start = 1.0;
      stop = 2.0;
      nsteps = 10;
      step = (stop-start)/float(nsteps);
      abs = false;
    } else {
      continue;
    }

    g_eff_MC.Set(nsteps);
    g_eff_Data.Set(nsteps);
    g_eff_residual.Set(nsteps);

    float min, max;
    int bin_inf, bin_sup;
    bin_inf = bin_sup = -1;
      
    int point=0;
    for(float cut=start; cut<=stop; cut+=step) {
      if(abs) {
        min = -1 * cut;
        max = cut;
      } else {
        min = start;
        max = cut;
      }
      bin_inf = histos_data[i]->FindBin(min);
      bin_sup = histos_data[i]->FindBin(max);
    
    
      float norm_MC, norm_data, sel_MC, sel_data, eff_MC, eff_data;
      
      norm_MC =  histos_MC[i]->Integral(1,histos_MC[i]->GetNbinsX());
      norm_data =  histos_data[i]->Integral(1, histos_data[i]->GetNbinsX());
      
      sel_MC =  histos_MC[i]->Integral(bin_inf, bin_sup);
      sel_data =  histos_data[i]->Integral(bin_inf, bin_sup);
      
      eff_MC = sel_MC / norm_MC;
      eff_data = sel_data / norm_data;
      
      float err_eff_MC = sqrt(eff_MC*(1-eff_MC)/fabs(histos_MC[i]->Integral()));
      float err_eff_data = sqrt(eff_data*(1-eff_data)/fabs(histos_data[i]->Integral()));

      g_eff_MC.SetPoint(point, cut, eff_MC);      
      g_eff_MC.SetPointError(point, 0., err_eff_MC);      

      g_eff_Data.SetPoint(point, cut, eff_data);      
      g_eff_Data.SetPointError(point, 0., err_eff_data);      
      
      g_eff_residual.SetPoint(point, cut, eff_data - eff_MC);      
      g_eff_residual.SetPointError(point, 0., sqrt(err_eff_data*err_eff_data + err_eff_MC*err_eff_MC));      

//       std::cout << "CUT = " << min << " - " << max << std::endl;
//       std::cout << "iecal = " << iecal << "\tvar = " << histos_data[i]->GetName() << "\t" <<  histos_MC[i]->GetName()
//                 << "\teff_MC = " << eff_MC << " +/- " << err_eff_MC
//                 << "\teff_data = " << eff_data << " +/- " << err_eff_data << std::endl;
//       std::cout << "Sel MC = " << sel_MC << "\tNorm MC = " << norm_MC 
//                 << "\tSel data = " << sel_data << "\t norm data = " << norm_data << std::endl;
      point++;
    } // loop over cuts

    // cosmetics
    g_eff_MC.SetLineColor(1);
    g_eff_MC.SetMarkerColor(1);
    g_eff_Data.SetLineColor(2);
    g_eff_Data.SetMarkerColor(2);

    TString XaxisTitle;
    if(name.Contains("dPhi")) XaxisTitle = "|#Delta #phi_{in}| cut";
    if(name.Contains("dEta")) XaxisTitle = "|#Delta #eta_{in}| cut";
    if(name.Contains("sigma")) XaxisTitle = "#sigma_{#eta#eta} cut";

    g_eff_Data.GetXaxis()->SetTitle(XaxisTitle);
    g_eff_Data.GetYaxis()->SetTitle("efficiency");
    
    g_eff_residual.GetXaxis()->SetTitle(XaxisTitle);
    g_eff_residual.GetYaxis()->SetTitle("Data - MC");

    TLegend* leg = new TLegend(0.15,0.75,0.40,0.90);
    leg->SetFillStyle(0); leg->SetBorderSize(0); leg->SetTextSize(0.03);
    leg->SetFillColor(0);
    leg->AddEntry(&g_eff_Data,"data","pl");
    leg->AddEntry(&g_eff_MC,"MC","pl");

    TCanvas c1;
    c1.Divide(1,2);
    c1.cd(1);
    g_eff_Data.Draw("APE");
    g_eff_MC.Draw("PE");
    leg->Draw();
    c1.cd(2);
    g_eff_residual.Draw("APE");
    TString namePS = name+"_eff.eps";
    c1.SaveAs(namePS);
  } // loop over variables


}



