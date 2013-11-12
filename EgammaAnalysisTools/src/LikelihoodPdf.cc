#include "EgammaAnalysisTools/include/LikelihoodPdf.h"
#include <stdlib.h>
#include <iostream>




LikelihoodPdf::LikelihoodPdf(const char* name, const char* species, int ecalsubdet, int ptbin) {
  _name = std::string(name);
  _species = std::string(species);
  _ecalsubdet = ecalsubdet;
  _ptbin = ptbin;
}



LikelihoodPdf::~LikelihoodPdf() {
}



void 
LikelihoodPdf::split(std::map<std::string,float> splitFractions, 
		     bool splitPdf) {

  char buffer[100];
  //! use different a-priori probabilities and different PDFs 
  //! depending by category
  if(splitFractions.size()>0 && splitPdf) {
    std::map<std::string,float>::const_iterator splitCatItr;
    for(splitCatItr=splitFractions.begin();splitCatItr!=splitFractions.end();splitCatItr++) {
      sprintf(buffer,"%sClass_%s_subdet%d_ptbin%d_%s",_name.c_str(),_species.c_str(),_ecalsubdet,_ptbin,splitCatItr->first.c_str());
      std::string totPdfName = std::string(buffer);
      _splitRule.insert( std::make_pair(splitCatItr->first,totPdfName) );
    }
  }

  //! use different a-priori, but same PDFs for all categories
  else if(splitFractions.size()>0) {
    std::map<std::string,float>::const_iterator splitCatItr;
    for(splitCatItr=splitFractions.begin();splitCatItr!=splitFractions.end();splitCatItr++) {
      sprintf(buffer,"%sUnsplit_%s_subdet%d_ptbin%d",_name.c_str(),_species.c_str(),_ecalsubdet,_ptbin);
      std::string totPdfName = std::string(buffer);
      _splitRule.insert( std::make_pair(splitCatItr->first,totPdfName) );
    }
  }

  //! do not split at all (same PDF's, same a-priori for all categories)
  else {
      sprintf(buffer,"%sUnsplit_%s_subdet%d_ptbin%d",_name.c_str(),_species.c_str(),_ecalsubdet,_ptbin);
      std::string totPdfName = std::string(buffer);
      _splitRule.insert( std::make_pair("NOSPLIT",totPdfName) );
  }
}


void 
LikelihoodPdf::initFromFile(TDirectory *dir) {
  _tdirectory = dir;
  std::map<std::string,std::string>::const_iterator ruleItr;
  for(ruleItr=_splitRule.begin();ruleItr!=_splitRule.end();ruleItr++) {
    TH1F *histo = (TH1F*)_tdirectory->Get(ruleItr->second.c_str());
    if(!histo) {
      std::cout << "LikelihoodPdf ERROR\thisto " << ruleItr->second.c_str() << " not found" << std::endl;
      exit(0);
    }
    else {
      _splitPdf.insert( std::make_pair(ruleItr->first,histo) );
    }
  }
}


float 
LikelihoodPdf::getVal(float x, std::string gsfClass, 
		      bool normalized) {
  _tdirectory->cd();
  TH1F *thePdf=0;
  if(_splitPdf.size()>1) thePdf=_splitPdf.find(gsfClass)->second;
  else thePdf=_splitPdf.find("NOSPLIT")->second;
  
  float prob=-1;

  if(normalized)
    prob=thePdf->GetBinContent(findBin(x,thePdf))/thePdf->Integral(0,thePdf->GetNbinsX()+1);
  else
    prob=thePdf->GetBinContent(findBin(x,thePdf));

  return prob;
}

int LikelihoodPdf::findBin(float x, TH1F* histo) {
  if (x<histo->GetXaxis()->GetXmin())
    return 0;
  else if (x>histo->GetXaxis()->GetXmax())
    return (histo->GetNbinsX()+1);
  else
    return histo->GetXaxis()->FindFixBin(x);
}
