#ifndef MONITOR_HH
#define MONITOR_HH

#include <vector>
#include <string>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TFile.h>

class Monitor {

public:

  Monitor(){};
  Monitor(int *nCand);
  Monitor(int *nCand, std::vector<int> *bestElements);
  virtual ~Monitor();

  // modifiers
  void setPath(std::string path);

  // methods

  // book 1 1D histo with variable var1. It should be a vector, so it will choose the best elements and the fake elements
  void book1D(std::string name, std::string title, float* var1, int nbins, float min, float max, std::string types="All", std::string opts="");
  // add a given variable to a reduced tree
  void Add2Tree(const char *name, float *var1);
  void DumpTree(const char *filename);

  void Fill(float weight=1.0);
  void FillTree();
  void ExcludeElements(std::vector<int> *excludeElements);
  void WritePs(std::string namePs);
  void WriteRoot(TFile *file);

private:
  std::string _path;
  
  // 1D histos
  std::vector<TH1F*> _All1D;
  std::vector<TH1F*> _Best1D;
  std::vector<TH1F*> _Fake1D;
  std::vector<std::string> _Opt;

  // Variables
  std::vector<float*> _Var1D;
  int *_nCand;
  std::vector<int> *_bestElements;
  std::vector<int> *_excludeElements;
  TTree *_data;
  
};

#endif
