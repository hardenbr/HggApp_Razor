#include "CommonTools/include/FiguresOfMeritEvaluator.h"

#include <iostream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TFile.h"
#include <math.h>
#include <stdio.h>

using namespace std;

FiguresOfMeritEvaluator::FiguresOfMeritEvaluator() {

  m_xmin = 0.0;
  m_xmax = 1.0;
  m_ymin = 0.0;
  m_ymax = 1.0;

  m_signalTitle = "signal";
  m_backgroundTitle = "background";

  _labelFont        = 42;
  _labelOffset      = 0.015;
  _axisLabelSize    = 0.050;
  _titleOffset      = 1.6;
  
}

void FiguresOfMeritEvaluator::addSignal(const char *nameVar, TH1F* sig) {
  m_signalHisto.push_back(sig);
  m_names.push_back(TString(nameVar));
}

void FiguresOfMeritEvaluator::addSignal(const char *nameVar, TH2F* sig) {
  m_signalHisto2D.push_back(sig);
  m_names2D.push_back(TString(nameVar));
}

void FiguresOfMeritEvaluator::addBackgrounds(TH1F* bkg0, TH1F* bkg1, 
					     TH1F* bkg2, TH1F* bkg3) {
  if(bkg1) bkg0->Add(bkg1);
  if(bkg2) bkg0->Add(bkg2);
  if(bkg3) bkg0->Add(bkg3);
  m_bkgHisto.push_back(bkg0);
}

void FiguresOfMeritEvaluator::addBackgrounds(TH2F* bkg0, TH2F* bkg1, 
					     TH2F* bkg2, TH2F* bkg3) {
  if(bkg1) bkg0->Add(bkg1);
  if(bkg2) bkg0->Add(bkg2);
  if(bkg3) bkg0->Add(bkg3);
  m_bkgHisto2D.push_back(bkg0);
}


TGraphErrors* FiguresOfMeritEvaluator::getFOM1D(const char *nameVar, int option) {
  
  TGraphErrors *outGraph = new TGraphErrors();

  int indexVar = -1;
  for(unsigned int ivar=0; ivar<m_signalHisto.size(); ivar++) {

    if (m_names[ivar].Contains(nameVar) && 
	TString(nameVar).Contains(m_names[ivar])) indexVar=ivar;
    
  }

  if ( indexVar==-1 ) {

    std::cout << "ERROR! The requested variable ( "
	      << nameVar << " ) is not in the list of known variables!" << std::endl;
    return 0;

  }

  TH1F *signal = m_signalHisto[indexVar];
  TH1F *background = m_bkgHisto[indexVar];
  const char *cutDir = m_direction[indexVar];

  if( signal && background ) {

    std::cout << "Integral of signal histogram " << signal->GetName() << " = " 
	      << signal->Integral() << std::endl;
    std::cout << "Integral of background histogram " << background->GetName() << " = " 
	      << background->Integral() << std::endl;
    
    TAxis *axisS = signal->GetXaxis();
    TAxis *axisB = background->GetXaxis();
    int nBinsSig = axisS->GetNbins();
    int nBinsBkg = axisB->GetNbins();
    
    if( nBinsSig!=nBinsBkg ) {
      std::cout << "ERROR! signal and background histograms have different binning." 
		<< std::endl;
      return 0;
    }

    // needed also overflow + underflow
    outGraph->Set(nBinsSig+2);

    double signalIntegral = signal->Integral(0,nBinsSig+1);
    double backgroundIntegral = background->Integral(0,nBinsSig+1);

    double tmpSignalIntegral=0.0;
    double tmpBackgroundIntegral=0.0;

    for ( int ibin=0; ibin<=nBinsSig+1; ibin++) {

      if( strcmp(cutDir,"<")==0 ) {
	tmpSignalIntegral = signal->Integral(0,ibin);
	tmpBackgroundIntegral = background->Integral(0,ibin);
      }
      else if( strcmp(cutDir,">")==0 ) {
	tmpSignalIntegral = signal->Integral(ibin,nBinsSig+1);
	tmpBackgroundIntegral = background->Integral(ibin,nBinsSig+1);
      } else if( strcmp(cutDir,"=")==0 ) {
	tmpSignalIntegral = signal->GetBinContent(ibin);
	tmpBackgroundIntegral = background->GetBinContent(ibin);
      } else {
	std::cout << "CONFIGURATION ERROR! direction of the cut not set." << std::endl
		  << "Please use: \">\" for var>x0 or  \"<\" for var<x0" << std::endl;
	return 0;
      }

      double signalEff = tmpSignalIntegral / signalIntegral;
      double backgroundEff = tmpBackgroundIntegral / backgroundIntegral;

      // std::cout << "bin =" << ibin << " sigEff=" << signalEff << " bkgEff=" << backgroundEff << std::endl;

      if( option == 0 ) {
	outGraph->SetPoint(ibin,signalEff,1-backgroundEff);
	outGraph->SetPointError(ibin,0,0);
      }
      else if( option == 1 ) {
	outGraph->SetPoint(ibin,backgroundEff,signalEff);
	double backgroundEffErr = sqrt(backgroundEffErr*(1-backgroundEffErr)/backgroundIntegral);
	outGraph->SetPointError(ibin,backgroundEffErr,0.);
      }
      else {
	std::cout << "unrecognized option" << std::endl;
	return 0;
      }

    }
  }
  
  else {
    std::cout << "ERROR! Cannot find signal or background histogram for variable "
	      << nameVar << std::endl;
    return 0;
  }

  return outGraph;

}

TGraphErrors* FiguresOfMeritEvaluator::getFOM2D(const char *nameVar, int option) {
  
  TGraphErrors *outGraph = new TGraphErrors();
  TGraphErrors *outGraphEnvelope = new TGraphErrors();

  int indexVar = -1;
  for(unsigned int ivar=0; ivar<m_signalHisto2D.size(); ivar++) {

    if (m_names2D[ivar].Contains(nameVar) && 
	TString(nameVar).Contains(m_names2D[ivar])) indexVar=ivar;
    
  }

  if ( indexVar==-1 ) {

    std::cout << "ERROR! The requested variable ( "
	      << nameVar << " ) is not in the list of known variables!" << std::endl;
    return 0;

  }

  TH2F *signal = m_signalHisto2D[indexVar];
  TH2F *background = m_bkgHisto2D[indexVar];
  TString cutDirXY = TString(m_directionXY[indexVar]);

  TObjArray* selectionTokens = cutDirXY.Tokenize(":");
  if (selectionTokens->GetEntries()!=2) {
    std::cout << "Wrong directions. Format should be like <:< " << selectionTokens->GetEntries() << std::endl;
    return 0;
  }
  TString cutDirX =((TObjString*)(*selectionTokens)[0])->GetString();
  TString cutDirY =((TObjString*)(*selectionTokens)[1])->GetString();

  if( signal && background ) {

    std::cout << "Integral of signal histogram " << signal->GetName() << " = " 
	      << signal->Integral() << std::endl;
    std::cout << "Integral of background histogram " << background->GetName() << " = " 
	      << background->Integral() << std::endl;
    
    TAxis *XaxisS = signal->GetXaxis();
    TAxis *YaxisS = signal->GetYaxis();
    TAxis *XaxisB = background->GetXaxis();
    TAxis *YaxisB = background->GetYaxis();
    int nXBinsSig = XaxisS->GetNbins();
    int nYBinsSig = YaxisS->GetNbins();
    int nXBinsBkg = XaxisB->GetNbins();
    int nYBinsBkg = YaxisB->GetNbins();
    
    if( nXBinsSig!=nXBinsBkg || nYBinsSig!=nYBinsBkg ) {
      std::cout << "ERROR! signal and background histograms have different binning." 
		<< std::endl;
      return 0;
    }

    // needed also overflow + underflow
    outGraph->Set((nXBinsSig+2)*(nYBinsSig+2));

    double signalIntegral = signal->Integral(0,nXBinsSig+1,0,nYBinsSig+1);
    double backgroundIntegral = background->Integral(0,nXBinsSig+1,0,nYBinsSig+1);

    double tmpSignalIntegral=0.0;
    double tmpBackgroundIntegral=0.0;

    int ibin=0;
    for ( int ixbin=0; ixbin<=nXBinsSig+1; ixbin++) {
      for ( int iybin=0; iybin<=nYBinsSig+1; iybin++) {
	
      if( strcmp(cutDirX,"<")==0 && strcmp(cutDirY,"<")==0 ) {
	tmpSignalIntegral = signal->Integral(0,ixbin,0,iybin);
	tmpBackgroundIntegral = background->Integral(0,ixbin,0,iybin);
      }
      else if( strcmp(cutDirX,"<")==0 && strcmp(cutDirY,">")==0 ) {
	tmpSignalIntegral = signal->Integral(0,ixbin,iybin,nYBinsSig+1);
	tmpBackgroundIntegral = background->Integral(0,ixbin,iybin,nYBinsSig+1);
      }
      else if( strcmp(cutDirX,">")==0 && strcmp(cutDirY,"<")==0 ) {
	tmpSignalIntegral = signal->Integral(ixbin,nXBinsSig+1,0,iybin);
	tmpBackgroundIntegral = background->Integral(ixbin,nXBinsSig+1,0,iybin);
      }
      else if( strcmp(cutDirX,">")==0 && strcmp(cutDirY,">")==0 ) {
	tmpSignalIntegral = signal->Integral(ixbin,nXBinsSig+1,iybin,nYBinsSig+1);
	tmpBackgroundIntegral = background->Integral(ixbin,nXBinsSig+1,iybin,nYBinsSig+1);
      }
      else {
	std::cout << "CONFIGURATION ERROR! direction of the cut not set." << std::endl
		  << "Please use: \">\" for var>x0 or  \"<\" for var<x0" << std::endl;
	return 0;
      }
      
      double signalEff = tmpSignalIntegral / signalIntegral;
      double backgroundEff = tmpBackgroundIntegral / backgroundIntegral;

      // std::cout << "bin =" << ibin << " sigEff=" << signalEff << " bkgEff=" << backgroundEff << std::endl;

      if( option == 0 ) {
	outGraph->SetPoint(ibin,signalEff,1-backgroundEff);
	outGraph->SetPointError(ibin,0,0);
      }
      else if( option == 1 ) {
	outGraph->SetPoint(ibin,signalEff,backgroundEff);
	double backgroundEffErr = sqrt(backgroundEffErr*(1-backgroundEffErr)/backgroundIntegral);
	outGraph->SetPointError(ibin,0.,backgroundEffErr);
      }
      else {
	std::cout << "unrecognized option" << std::endl;
	return 0;
      }
      
      ibin++;
      }
    }

    // now get the envelope of the 2D
    Double_t* x=outGraph->GetX();
    Double_t* y=outGraph->GetY();
    int npoints=outGraph->GetN();
    
    //float effstep=1.0/float(nXBinsSig);
    float effstep = 0.03;
    int ebin=0;
    for(float eff=m_xmin;eff<m_xmax;eff+=effstep) {
      // look for the min bkg efficiency compatible with the efficiency of the step
      float ymin=1.0;
      for(int i=0;i<npoints;i++) {
	if(x[i]>eff-effstep && x[i]<=eff && y[i]<ymin) ymin=y[i]; 
      }
      outGraphEnvelope->SetPoint(ebin,eff,ymin);
      ebin++;
    }

  } else {
    std::cout << "ERROR! Cannot find signal or background histogram for variable "
	      << nameVar << std::endl;
    return 0;
  }
  
  return outGraphEnvelope;

}


void FiguresOfMeritEvaluator:: drawResults(const char *fileName, int option) {

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  if( m_signalHisto.size()!=m_bkgHisto.size() ) {
    std::cout << "ERROR! for some variable or signal or background histo is missing. Exiting!" << std::endl;
    return;
  }
  
  TPaveText *text = new TPaveText(0.15,0.90,0.77,0.98,"brNDC");
  text->AddText("CMS Preliminary         #sqrt{s} = 8 TeV,  L = 19.6 fb^{-1}");
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(132);
  text->SetTextSize(0.04);

  TCanvas c1("c1","",600,600);
  c1.Range(-1.146789,-2319.078,5.688073,12419.95);
  c1.SetFillColor(0);
  c1.SetBorderMode(0);
  c1.SetBorderSize(2);
  c1.SetLeftMargin(0.1677852);
  c1.SetFrameBorderMode(0);
  c1.SetFrameBorderMode(0);
  c1.cd();

  float legxmin, legxmax, legymin,legymax;
  if(option==0) {
    legxmin=0.20;
    legxmax=0.40;
    legymin=0.20;
    legymax=0.30;
  } else {
    legxmin=0.50;
    legxmax=0.70;
    legymin=0.20;
    legymax=0.40;
  }

  TLegend* leg = new TLegend(legxmin,legymin,legxmax,legymax);
  leg->SetBorderSize(     0);
  leg->SetFillColor (     0);
  leg->SetTextAlign (    12);
  leg->SetTextFont  (_labelFont);
  leg->SetTextSize  (  0.03);
    
  std::vector<TGraphErrors*> graphs;

  // loop over 1D
  for( unsigned int ivar=0; ivar<m_signalHisto.size(); ivar++) {
    const char *name = m_names[ivar].Data();
    std::cout << "---> processing " << name << "..."; 
    TGraphErrors *graph = getFOM1D(name,option);
    graphs.push_back(graph);
    leg->AddEntry(graph,name,"p");
  }

  // loop over 2D
  for( unsigned int ivar=0; ivar<m_signalHisto2D.size(); ivar++) {
    const char *name = m_names2D[ivar].Data();
    std::cout << "---> processing " << name << "..."; 
    TGraphErrors *graph = getFOM2D(name,option);
    graphs.push_back(graph);
    leg->AddEntry(graph,name,"p");
  }

  // draw the results
  TString fullname(fileName);
  TObjArray *tokens = fullname.Tokenize(".");
  const char *basename = (((TObjString*)(*tokens)[0])->GetString()).Data();

  TString pdf = TString(basename)+TString(".pdf");
  TString png = TString(basename)+TString(".png");
  TString macro = TString(basename)+TString(".C");

  TString root = TString(basename)+TString(".root");
  TFile *tfile = TFile::Open(root,"recreate");

  for(int ig=0;ig<(int)graphs.size();++ig) {
    TGraphErrors* graph = graphs[ig];
    if( graph ) {

      char nameg[50];
      sprintf(nameg,"graph_%d",ig);
      graph->SetName(nameg);
      graph->SetTitle("");
      if(ig==0) graph->SetMarkerStyle(kFullDotMedium);
      else if(ig==1) graph->SetMarkerStyle(kFullSquare);
      else if(ig==2) graph->SetMarkerStyle(kFullTriangleDown);
      else if(ig==3) graph->SetMarkerStyle(kFullTriangleUp);
      else graph->SetMarkerStyle(kFullStar+ig);
      int defColor;
      if(ig==0) defColor=kRed+1;
      else if(ig==1) defColor=kAzure-6;
      else if(ig==2) defColor=kTeal+3;
      else if(ig==3) defColor=kViolet+3;
      else defColor = ig+1;
      // skip yellow
      if(defColor==5) defColor=kOrange+4;
      graph->SetMarkerColor(defColor);
      graph->SetMarkerSize(1.4);
      graph->SetLineColor(defColor);
      graph->SetLineWidth(2);
      graph->GetXaxis()->SetRangeUser(m_xmin,m_xmax);
      graph->GetYaxis()->SetRangeUser(m_ymin,m_ymax);
      
      std::string sigSuffix = " efficiency";
      std::string xAxisName = std::string(m_backgroundTitle) + sigSuffix;
      
      std::string bkgSuffix = (option==0) ? " rejection" : " efficiency";
      std::string yAxisName = std::string(m_signalTitle) + bkgSuffix;
      
      graph->GetXaxis()->SetTitle(xAxisName.c_str());
      graph->GetYaxis()->SetTitle(yAxisName.c_str());
      
      AxisFonts(graph->GetXaxis(), "x", graph->GetXaxis()->GetTitle());
      AxisFonts(graph->GetYaxis(), "y", graph->GetYaxis()->GetTitle());
      
      if(ig==0) graph->Draw("APE2");
      else  graph->Draw("PE2");

      tfile->cd();
      graph->Write();

    }

  }
  
  leg->Draw();
  text->Draw();
  tfile->Close();

  c1.SaveAs(pdf);
  c1.SaveAs(png);
  c1.SaveAs(macro);
  
}

void FiguresOfMeritEvaluator::setRange(double xmin, double xmax, double ymin, double ymax) {

  m_xmin = xmin;
  m_xmax = xmax;
  m_ymin = ymin;
  m_ymax = ymax;

}

  //------------------------------------------------------------------------------
  // AxisFonts
  //------------------------------------------------------------------------------
void FiguresOfMeritEvaluator::AxisFonts(TAxis*  axis,
					TString coordinate,
					TString title)
{
  axis->SetLabelFont  (_labelFont  );
  // axis->SetLabelOffset(_labelOffset);
  // axis->SetLabelSize  (_axisLabelSize);
  // axis->SetNdivisions (  505);
  axis->SetTitleFont  (_labelFont);
  axis->SetTitleOffset(  1.5);
  // axis->SetTitleSize  (_axisLabelSize);
  if (coordinate == "y") axis->SetTitleOffset(1.7);
  if (coordinate == "x") axis->SetTitleOffset(1.2);
  
  axis->SetTitle(title);
}
