#include <string>
#include <iostream>

#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TPostScript.h>

#include "CommonTools/include/Monitor.hh"

Monitor::Monitor(int *nCand) {
  _nCand=nCand;

  // an empty default vector
  std::vector<int> *defBestElements=new std::vector<int>;
  _bestElements=defBestElements;

  std::vector<int> *defExcludeElements=new std::vector<int>;
  _excludeElements=defExcludeElements;

  _data = new TTree("data","data");
}

Monitor::Monitor(int *nCand, std::vector<int> *bestElements) {
  _nCand=nCand;
  _bestElements=bestElements;
}

Monitor::~Monitor() {}

void Monitor::setPath(std::string path) {
  _path=path;
}

void Monitor::book1D(std::string name, std::string title, float* var1, int nbins, float min, float max, std::string types, std::string opts) {
  TH1F* histo = new TH1F("histo",title.c_str(),nbins,min,max);
  
  std::string::size_type loc=types.find("Best");
  if(loc!=std::string::npos) {
    _Best1D.push_back((TH1F*)histo->Clone());
    _Best1D.back()->SetName((name+"_Best").c_str());
    _Best1D.back()->SetTitle("");
    _Best1D.back()->GetXaxis()->SetTitle(title.c_str());
  }

  loc=types.find("Fake");
  if(loc!=std::string::npos) {
    _Fake1D.push_back((TH1F*)histo->Clone());
    _Fake1D.back()->SetName((name+"_Fake").c_str());
    _Fake1D.back()->SetTitle("");
    _Fake1D.back()->GetXaxis()->SetTitle(title.c_str());
  }

  loc=types.find("All");
  if(loc!=std::string::npos) {
    _All1D.push_back((TH1F*)histo->Clone());
    _All1D.back()->SetName((name+"_All").c_str());
    _All1D.back()->SetTitle("");
    _All1D.back()->GetXaxis()->SetTitle(title.c_str());
  }

  _Var1D.push_back(var1);
  _Opt.push_back(opts);

}

void Monitor::Add2Tree(const char* name, float *var1) {
  std::string m(name);
  const char *fullname=((m+"/").c_str());
  _data->Branch(name,var1,fullname);
}

void Monitor::DumpTree(const char *filename) {
  TFile rootfile(filename,"RECREATE");
  _data->Write();
  rootfile.Close();
}

void Monitor::ExcludeElements(std::vector<int> *excludeElements) {
  _excludeElements=excludeElements;
}

void Monitor::Fill(float weight) {
  // fill the histograms
  for(int i=0; i<(int)_All1D.size(); i++) {
    float* var1 = _Var1D[i];
    if(_nCand!=0) {
      for(int j=0; j<*_nCand; j++) {
	_All1D[i]->Fill(var1[j],weight);
	// look for the best candidates index
	bool isOneOfBest=false;
	std::vector<int>::const_iterator iterBests;
	for(iterBests=_bestElements->begin(); iterBests!=_bestElements->end(); iterBests++) {
	  if(j==*iterBests) { 
	    isOneOfBest=true; 
	    break;
	  }
	}

	// look for candidates to exclude index
 	bool isOneOfExcluded=false;
 	if(_excludeElements) {
 	  std::vector<int>::const_iterator iterExcluded;
	  for(iterExcluded=_excludeElements->begin(); iterExcluded!=_excludeElements->end(); iterExcluded++) {
	    if(j==*iterExcluded) { 
	      isOneOfExcluded=true; 
	      break;
	    }
	  }
 	}
	if(isOneOfBest && _Best1D.size()!=0 && _bestElements->size()!=0 && !isOneOfExcluded )  _Best1D[i]->Fill(var1[j],weight);
	if(!isOneOfBest && _Fake1D.size()!=0 && _bestElements->size()!=0 && !isOneOfExcluded ) _Fake1D[i]->Fill(var1[j],weight);
      }
    }
    else {
      _All1D[i]->Fill(var1[0],weight);
    }
  }
}

void Monitor::FillTree() {
  // fill the dataset
  _data->Fill();
}

void Monitor::WritePs(std::string namePs) {
  int npad=0;
  bool best(false),fake(false),all(false);
  if(_Best1D.size()!=0) {best=true; npad++;}
  if(_Fake1D.size()!=0) {fake=true; npad++;}
  if(_All1D.size()!=0)  {all=true; npad++;}

  ///////////// Default definitions /////////////////////////
  gROOT->SetBatch(kTRUE);
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  //  gStyle->SetOptStat(1111111);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptStat(0);  // Show nothing
  gStyle->SetStatStyle(0001); // for a completely transparent stat box
  gStyle->SetOptFit(111110); 
  gStyle->SetStatX(0.85);
  gStyle->SetStatY(0.95);
  gStyle->SetOptFile(1); 
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(.1);
  gStyle->SetMarkerColor(1);
  gStyle->SetTitleBorderSize(0);  // no border around histogram title (font size can't be calculated anyways ...)
  TCanvas c0("c0","--c0--",356,0,656,800);
  c0.ToggleEventStatus();
  TLatex tl;
  tl.SetNDC(kTRUE);
  TPostScript *ps = new TPostScript(namePs.c_str() ,111);
  ///////////////////////////////////////////////////////////

  if(npad==2) c0.Divide(1,2);  
  if(npad==3) c0.Divide(1,3);  
  if(npad==4) c0.Divide(2,2);


  //Draw 1D
  for(int i=0;i<int(_All1D.size());i++) {
    ps->NewPage();
    bool log(false);
    if(_Opt[i]=="Log") log=true;
    int j=0;
    c0.cd(++j);
    if(_Best1D.size()!=0){
      gPad->SetLogy(log);
      _Best1D[i]->Draw();
    }
    c0.cd(++j);
    if(_Fake1D.size()!=0){
      gPad->SetLogy(log);
      _Fake1D[i]->Draw();
    }
    c0.cd(++j);
    if(_All1D.size()!=0){
      gPad->SetLogy(log);
      _All1D[i]->Draw();
    }
    //    c0.cd(++j);
    c0.Update();
  }
  ps->Close();
  std::cout << "*** Saved histograms in: " << namePs << " ***" << std::endl;
}

void Monitor::WriteRoot(TFile *file) {
  file->mkdir(_path.c_str());
  file->cd(_path.c_str());
  for(int i=0; i<int(_All1D.size());i++){
    if(_Best1D.size()!=0){
      _Best1D[i]->SetName(_Best1D[i]->GetName());
      _Best1D[i]->Write();
    }
    if(_Fake1D.size()!=0){
      _Fake1D[i]->SetName(_Fake1D[i]->GetName());
      _Fake1D[i]->Write();
    }
    if(_All1D.size()!=0){
      _All1D[i]->SetName(_All1D[i]->GetName());
      _All1D[i]->Write();
    }
  }
}
