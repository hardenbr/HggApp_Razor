#include "CommonTools/include/Counters.hh"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>

void Counters::AddVar(string name) {
  _names.push_back(name); 
  _counts.push_back(int(0));
}

void Counters::IncrVar(string name) {
  for(int i=0;i<int(_names.size());i++){
    if(_names[i]==name) _counts[i]++;
  }
}

void Counters::IncrVar(string name, float weight) {
  for(int i=0;i<int(_names.size());i++){
    if(_names[i]==name) _counts[i]+=weight;
  }
}

float Counters::GetVar(string name) {
  for(int i=0;i<int(_names.size());i++){
    if(_names[i]==name) return (float)_counts[i];
  }
}

void Counters::SetTitle(const char* title) {
  _title=std::string(title);
}

void Counters::IncrNevents(){
  _nevents++;
}

void Counters::Draw(){
  cout<<"*****************************************"<<endl;
  cout<<"*------- Counters output table ----------"<<endl;
  cout<<"*-------  for " << _title << " ----------"<<endl;
  cout<<"*-------       Efficiencies    ----------"<<endl;
  cout<<"*"<<endl;
  std::cout<<std::setiosflags( std::ios::fixed );     // Format x.y 
  std::cout<<std::setiosflags( std::ios::showpoint ); // 0.10 
  std::cout<<std::setprecision(3); //2 decimal places

  for(int i=0;i<int(_names.size());i++){
    //float eff=float(_counts[i])/_nevents;
    //cout<<"* "<<_names[i]<<":\t"<<eff
    //	<<" +/- "<<sqrt(eff*(1-eff)/_nevents)<<endl;

    cout<<"* "<<_names[i]<<":\t"<<_counts[i]<< endl;
  }
  cout<<"*"<<endl;
  cout<<"*****************************************"<<endl;
}

void Counters::Draw(string name1,string name2) {
  
  std::cout<<std::setiosflags( std::ios::fixed );     // Format x.y 
  std::cout<<std::setiosflags( std::ios::showpoint ); // 0.10 
  std::cout<<std::setprecision(3); //2 decimal places
  
  for(int i=0;i<int(_names.size());i++) {
    for (int j=0;j<int(_names.size());j++) {
      if(_names[i]==name1 && _names[j]==name2) {
	float eff=float(_counts[i])/float(_counts[j]);
	cout<<"-*- "<<_names[i]<<":\t"<< 100 * eff
	    <<" +/- "<< 100 * sqrt(eff*(1-eff)/float(_counts[j]))<<endl;
      }
    }
  }
}

void Counters::Save(const char* filename, const char* option) {
  
  TFile *file = TFile::Open(filename,option);
  file->cd();

  TString branchTitle(_title.c_str());
  branchTitle.ReplaceAll(" ","_");
  branchTitle.ReplaceAll("+","_");

  TTree tree(branchTitle,"Efficiency Counters");

  int ncuts = _counts.size();
  float counts[ncuts];

  tree.Branch("nCuts",&ncuts,"nCuts/I");
  tree.Branch("nSel",counts,"nSel[nCuts]/F");

  for(int i=0;i<ncuts;i++) counts[i]=_counts[i];

  tree.Fill();
  tree.Write();
  file->Close();

}
