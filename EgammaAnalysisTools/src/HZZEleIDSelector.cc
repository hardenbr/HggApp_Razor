#include "EgammaAnalysisTools/include/HZZEleIDSelector.hh"
#include "EgammaAnalysisTools/include/eIDMVACuts.h"

#include <math.h>
#include <iostream>

using namespace std;

int HZZEleIDSelector::etabin(float eta) {
  if(fabs(eta)<0.8) return 0;
  else if(fabs(eta)<1.479) return 1;
  else return 2;
}

int HZZEleIDSelector::ptbinTrg(float pt) {
  if(pt<20) return 0;
  else return 1;
}

int HZZEleIDSelector::ptbinNoTrg(float pt) {
  if(pt<10) return 0;
  else return 1;
}

bool HZZEleIDSelector::output(float pt, float eta, float bdt, float iso, 
			      HZZEleIDSelector::wpfulliso WP, HZZEleIDSelector::mvatype type,
                              HZZEleIDSelector::cutblock cuts) {
  int etab=etabin(eta);
  int ptb=-1;
  if(type==kMVABiased) ptb=ptbinTrg(pt);
  else ptb=ptbinNoTrg(pt);
  float bdtcut=999.;
  float isocut=-999.;

  // special WP that has the same eff as HWW 2011
  if(WP==kWPHWW && type==kMVABiased) {
    bdtcut=cutbdtfulliso_hww[ptb][etab];
    isocut=cutfulliso_hww[ptb][etab];
    if(cuts==idonly) return (bdt>bdtcut);
    if(cuts==isoonly) return (iso/pt<isocut);
    return (bdt>bdtcut && iso/pt<isocut);
  }

  // special WP that has the same fake rate as HZZ 2011
  if((WP==kMVALoose || WP==kMVATight) && type==kMVAUnbiased) {
    int wp = (WP==kMVALoose) ? 0 : 1;
    bdtcut=cutbdtfulliso_hzz[ptb][etab][wp];
    isocut=cutmvaiso_hzz[ptb][etab][wp];
    if(cuts==idonly) return (bdt>bdtcut);  
    if(cuts==isoonly) return (iso>isocut);  // attention: this is MVA output, not a relative iso!
    return (bdt>bdtcut && iso>isocut);
  }

  // optimized working points
  if(type==kMVABiased) {
    bdtcut=cutbdtfulliso_biased[ptb][etab][WP];
    isocut=cutfulliso_biased[ptb][etab][WP];
  } else if(type==kMVAUnbiased) {
    bdtcut=cutbdtfulliso_unbiased[ptb][etab][WP];
    isocut=cutfulliso_unbiased[ptb][etab][WP];    
  } else {
    cout << " WARNING! type of MVA unset! " << endl;
  }
  // cout << "   eta = " << eta << " pt = " << pt 
  //      << "   bdt = " << bdt << " (cut is: " << bdtcut << ")" 
  //      << "   iso = " << iso << " (cut is: " << isocut << ")"
  //      << endl;
  if(cuts==idonly) return (bdt>bdtcut);  
  if(cuts==isoonly) return (iso/pt<isocut);
  return (bdt>bdtcut && iso/pt<isocut);
}

bool HZZEleIDSelector::output(float pt, float eta, float bdt, float iso, 
			      HZZEleIDSelector::wpchiso WP, HZZEleIDSelector::mvatype type,
			      HZZEleIDSelector::cutblock cuts) {
  int etab=etabin(eta);
  int ptb=-1;
  if(type==kMVABiased) ptb=ptbinTrg(pt);
  else ptb=ptbinNoTrg(pt);
  float bdtcut=999.;
  float isocut=-999.;
  if(type==kMVABiased) {
    bdtcut=cutbdtchiso_biased[ptb][etab][WP];
    isocut=cutchiso_biased[ptb][etab][WP];
  } else if(type==kMVAUnbiased) {
    bdtcut=cutbdtchiso_unbiased[ptb][etab][WP];
    isocut=cutchiso_unbiased[ptb][etab][WP];    
  } else {
    cout << " WARNING! type of MVA unset! " << endl;
  }
  if(cuts==idonly) return (bdt>bdtcut);
  if(cuts==isoonly) return (iso/pt<isocut);
  return (bdt>bdtcut && iso/pt<isocut);
}
