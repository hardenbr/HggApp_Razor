#include <math.h>

#include "EgammaAnalysisTools/include/eIDCiChzzCuts.h"
#include "EgammaAnalysisTools/include/eIDCiChzzSelector.hh"

bool eIDCiChzzSelector::compute_eid_cut(float x, float et, float cut_min, float cut_max, bool gtn) {

  float et_min = 10;
  float et_max = 40;

  bool accept = false;
  float cut = cut_max;    //  the cut at et=40 GeV

  if(et < et_max) {
    //cut = cut_max+(1/et-1/et_min)*(cut_min-cut_max)/(1/et_max-1/et_min);
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  }

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  //std::cout << x << " " << cut_min << " " << cut << " " << cut_max << " " << et << " " << accept << std::endl;
  return accept;
}

int eIDCiChzzSelector::ElectronClassification(Float_t myeta, Float_t fbrem, Float_t eopin, Int_t isPflow, bool specialCategories) {

  int cat = -1;
  float eta = fabs(myeta);

  if (eta < 1.479) {
    if (fbrem>=0.12 && eopin >0.9 && eopin < 1.2)                           // bremming barrel
      cat = 0;
    else if (specialCategories && ((eta >  .445   && eta <  .45  ) ||       // barrel crack electrons
				   (eta >  .79    && eta <  .81  ) ||
				   (eta > 1.137   && eta < 1.157 ) ||
				   (eta > 1.47285 && eta < 1.4744)))
      cat = 6;                                       
    else if (isPflow && specialCategories)                // tracker driven electrons
      cat = 8;
    else if (fbrem < 0.12)                                                  // low brem barrel 
      cat = 1;
    else                                             // bad track in barrel
      cat = 2;
  } else {
    if (fbrem>=0.2 && eopin >0.82 && eopin < 1.22)   // bremming endcap
      cat = 3;
    else if (specialCategories && (eta > 1.5 && eta < 1.58))                // endcap crack electrons
      cat = 7;                                       
    else if (specialCategories && isPflow)              // tracker driven electrons
      cat = 8;
    else if (fbrem < 0.2)                            // low brem endcap 
      cat = 4;
    else                                             // bad track in barrel
      cat = 5;
  }

  return cat;
}

Int_t eIDCiChzzSelector::ElectronId_V03(Float_t et, Float_t eta, Float_t sieie, Float_t eop, Float_t eseedopin, Float_t fbrem, 
					Float_t tkiso03, Float_t ecaliso04, Float_t hcaliso04, Float_t ip, 
					Float_t mishits, Float_t detain, Float_t dphiin, Float_t hoe, Float_t dcot,
					Float_t dist, Int_t isPflow, Int_t level, bool specialCategories) {
  
  int result = 0;

  int cat = ElectronClassification(eta, fbrem, eop, isPflow, specialCategories);
  const int ncuts = 10;
  bool cut_results[ncuts];

  float eseedopincor = eseedopin + fbrem;
  if(fbrem < 0)
    eseedopincor = eseedopin;

  float iso_sum = tkiso03 + ecaliso04 + hcaliso04;

  if(fabs(eta) > 1.5) 
    iso_sum += (fabs(eta) - 1.5)*1.09;

  float iso_sumoet = iso_sum*(40./et);

  Float_t dcotdist = ((0.04 - std::max(fabs(dist), fabs(dcot))) > 0?(0.04 - std::max(fabs(dist), fabs(dcot))):0);

  for (int cut=0; cut<ncuts; cut++) {
    float cut_min, cut_max;
    switch (cut) {
    case 0:
      cut_results[cut] = compute_eid_cut(fabs(detain), et, cutdetainlv03[0][level][cat], cutdetainv03[0][level][cat]);
      break;
    case 1:
      cut_results[cut] = compute_eid_cut(fabs(dphiin), et, cutdphiinlv03[0][level][cat], cutdphiinv03[0][level][cat]);
      break;
    case 2:
      cut_results[cut] = (eseedopincor > cuteseedopcorv03[0][level][cat]);
      break;
    case 3:
      cut_results[cut] = compute_eid_cut(hoe, et, cuthoelv03[0][level][cat], cuthoev03[0][level][cat]);
      break;
    case 4:
      cut_results[cut] = compute_eid_cut(sieie, et, cutseelv03[0][level][cat], cutseev03[0][level][cat]);
      break;
    case 5:
      cut_results[cut] = compute_eid_cut(iso_sumoet, et, cutiso_sumoetlv03[0][level][cat], cutiso_sumoetv03[0][level][cat]);
      break;
    case 6:
      cut_results[cut] = (iso_sum < cutiso_sumv03[0][level][cat]);
      break;
    case 7:
      cut_results[cut] = compute_eid_cut(fabs(ip), et, cutip_gsflv03[0][level][cat], cutip_gsfv03[0][level][cat]);
      break;
    case 8:
      cut_results[cut] = (mishits < cutfmishitsv03[0][level][cat]);
      break;
    case 9:
      cut_results[cut] = (dcotdist < cutdcotdistv03[0][level][cat]);
      break;
    }
  }
    
  // ID part
  if (cut_results[0] & cut_results[1] & cut_results[2] & cut_results[3] & cut_results[4])
    result = result + 1;
  
  // ISO part
  if (cut_results[5] & cut_results[6])
    result = result + 2;
  
  // IP part
  if (cut_results[7])
    result = result + 8;
  
  // Conversion part
  if (cut_results[8] & cut_results[9])
    result = result + 4;

  return result;
}

//cat = ElectronClassification(fabs(eta), fbrem, eop, isPflow, true);
//eId[level] = ElectronId_V03(et, fabs(eta), sieie, eop, eseedop, fbrem, tkiso, ecaliso, hcaliso, ip, mishits, detain, dphiin, hoe, dcot, dist, isPflow, level, true);
