#include "include/eIDSimpleCutsSelector.hh"
#include <math.h>

eIDSimpleCutsSelector::eIDSimpleCutsSelector(float deta, float dphi, float see, float hoe, 
					     float d0, float dz, float IoEmIoP) {
  deta_=deta;
  dphi_=dphi;
  see_=see;
  hoe_=hoe;
  d0_=d0;
  dz_=dz;
  IoEmIoP_=IoEmIoP;
}

bool eIDSimpleCutsSelector::output(float eta, int WP) {
  float cut_deta[2][4] = 
    {
      { 0.007, 0.007, 0.004, 0.004 }, // barrel
      { 0.01, 0.009, 0.007, 0.005}, // endcap
    };
  float cut_dphi[2][4] = 
    {
      { 0.8, 0.15, 0.06, 0.03 },
      {0.7, 0.10, 0.03, 0.03 }
    };
  float cut_see[2][4] = 
    {
      { 0.01, 0.01, 0.01, 0.01 },
      { 0.03, 0.03, 0.03, 0.02 }
    };
  float cut_hoe[2][4] =
    {
      { 0.15, 0.12, 0.12, 0.12 },
      { 1000, 0.10, 0.10, 0.10 }
    };
  float cut_d0[2][4] = 
    {
      { 0.04, 0.02, 0.02, 0.02 },
      { 0.04, 0.02, 0.02, 0.02 }
    };
  float cut_dz[2][4] =
    {
      { 0.2, 0.2, 0.1, 0.1 },
      { 0.2, 0.2, 0.1, 0.1 }
    };
  float cut_IoEmIoP[2][4] =
    {
      { 1000, 0.05, 0.05, 0.05 },
      { 1000, 0.05, 0.05, 0.05 }
    };
  
  int etabin = (fabs(eta)<1.479) ? 0 : 1;
  return fabs(deta_)<cut_deta[etabin][WP] &&
    fabs(dphi_)<cut_dphi[etabin][WP] &&
    see_<cut_see[etabin][WP] &&
    hoe_<cut_hoe[etabin][WP] &&
    fabs(d0_)<cut_d0[etabin][WP] &&
    fabs(dz_)<cut_dz[etabin][WP] &&
    fabs(IoEmIoP_)<cut_IoEmIoP[etabin][WP];

}
