#include "EgammaAnalysisTools/include/SimpleCutsIDSelector.hh"
#include "EgammaAnalysisTools/include/eIDSimpleCuts.h"

#include <math.h>
#include <iostream>

using namespace std;

int SimpleCutsIDSelector::etabin(float eta) {
  if(fabs(eta)<1.479) return 0;
  else return 1;
}

bool SimpleCutsIDSelector::outputid(float eta, float see, float dphi, float deta, float hoe, 
				    SimpleCutsIDSelector::wp WP) {
  int etab=etabin(eta);
  float seecut=cutid[etab][WP][simplecuts::seecut];
  float dphicut=cutid[etab][WP][simplecuts::dphicut];
  float detacut=cutid[etab][WP][simplecuts::detacut];
  float hoecut=cutid[etab][WP][simplecuts::hoecut];
  return (see<seecut && fabs(dphi)<dphicut && fabs(deta)<detacut && hoe<hoecut);
}

bool SimpleCutsIDSelector::outputiso(float eta, float tk, float ecal, float hcal,
				     SimpleCutsIDSelector::wp WP) {
  int etab=etabin(eta);
  float tkisocut=cutiso[etab][WP][simplecuts::tkcut];
  float ecalisocut=cutiso[etab][WP][simplecuts::ecalcut];
  float hcalisocut=cutiso[etab][WP][simplecuts::hcalcut];
  return (fabs(tk)<tkisocut && fabs(ecal)<ecalisocut && fabs(hcal)<hcalisocut);
}

bool SimpleCutsIDSelector::outputconv(float eta, int misshits, float dist, float dcot, 
				      SimpleCutsIDSelector::wp WP) {
  int etab=etabin(eta);
  float misshitscut=cutconvrej[etab][WP][simplecuts::misshitscut];
  float distcut=cutconvrej[etab][WP][simplecuts::distcut];
  float dcotcut=cutconvrej[etab][WP][simplecuts::dcotcut];
  return(misshits<misshitscut && (dist>distcut || dcot>dcotcut));
}
