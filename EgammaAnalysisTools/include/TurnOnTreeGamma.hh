#ifndef TurnOnTreeGamma_H
#define TurnOnTreeGamma_H

#include <TFile.h>
#include <TTree.h>

class TurnOnTreeGamma {
public:
  
  TurnOnTreeGamma(const char *filename);
  ~TurnOnTreeGamma();

  //! add more variables
  void addTandPinfo();
  void addHLTmatchInfo();
  void addMore();
  void addRunInfos();

  //! fill the tree with the wanted set of variables
  void fillVariables(float ptg, float etag, float phig, int idg, int hltmatchg);
  void fillTandPinfo(float mass, float pt1, float pt2, float eta1, float eta2, int tag1, int tag2, int isAZ, int isAZt);
  void fillHLTmatchInfo(int hlt20, int hlt30, int hlt50, int hlt75, int hlt90);
  void fillMore(float nVtx, float rho);
  void fillRunInfos(int run, int lumi, int event);

  void store();
  void save();

private:

  float myPtGamma, myEtaGamma, myPhiGamma;
  int myIdGamma;
  int myHltmatchGamma;

  float myZmass;
  float myPtEle1, myPtEle2;
  float myEtaEle1, myEtaEle2;
  int myTagEle1, myTagEle2;
  int myIsAZ, myIsAZ_tight; 

  int myPassHLT20, myPassHLT30, myPassHLT50, myPassHLT75, myPassHLT90; 
  
  int myRun, myLS, myEvent;
  float myNVtx, myRho;

  TFile *myFile;
  TTree *myTreeGamma;

};

#endif // TurnOnTreeGamma_H
