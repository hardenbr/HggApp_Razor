#ifndef FakeTree_H
#define FakeTree_H

#include <TFile.h>
#include <TTree.h>

class FakeTree {
public:
  
  FakeTree(const char *filename);
  ~FakeTree();
  
  void fill(float met, float masst, float me, float jpt, float dpt, float dphi);

  void store();
  void save();

private:

  float myMet, myWmt, myMee, myJetPt, myDenomPt, myDeltaPhi;
  
  TFile *myFile;
  TTree *myTree;

};

#endif // FakeTree_H
