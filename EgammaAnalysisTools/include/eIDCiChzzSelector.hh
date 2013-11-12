#ifndef HZZCICSELECTOR_H
#define HZZCICSELECTOR_H

#include "TROOT.h"

class eIDCiChzzSelector {
public:
  eIDCiChzzSelector() {};
  ~eIDCiChzzSelector() {};
  
  Int_t ElectronId_V03(Float_t et, Float_t eta, Float_t sieie, Float_t eop, Float_t eseedopin, Float_t fbrem, 
		       Float_t tkiso03, Float_t ecaliso04, Float_t hcaliso04, Float_t ip, 
		       Float_t mishits, Float_t detain, Float_t dphiin, Float_t hoe, Float_t dcot,
		       Float_t dist, Int_t isPflow, Int_t level, bool specialCategories);
  

private:

  bool compute_eid_cut(float x, float et, float cut_min, float cut_max, bool gtn=false);
  int ElectronClassification(Float_t myeta, Float_t fbrem, Float_t eopin, Int_t isPflow, bool specialCategories);

};

#endif

