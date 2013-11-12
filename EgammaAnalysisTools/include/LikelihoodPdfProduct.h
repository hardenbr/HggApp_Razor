#ifndef LikelihoodPdfProduct_h
#define LikelihoodPdfProduct_h

#include "EgammaAnalysisTools/include/LikelihoodSpecies.h"
#include "EgammaAnalysisTools/include/LikelihoodPdf.h"
#include <TDirectory.h>
#include <string>
#include <vector>
#include <map>

class LikelihoodPdfProduct {
 public:
  LikelihoodPdfProduct(const char* name, int ecalsubdet, int ptbin);
  ~LikelihoodPdfProduct();
  
  //! initialize the Pdf's from file
  void initFromFile(TDirectory *dir);

  //! add a species (hypothesis) to the likelihood, with a priori probability 
  void addSpecies(const char* name, float priorWeight=1.);

  //! add a pdf for a species, splitted or not
  void addPdf(const char* specname, const char* name, bool splitPdf=false);

  //! set the fraction of one category for a given species
  void setSplitFrac(const char* specname, const char* catName, float frac);

  //! get the likelihood ratio p(a priori) * L(specName) / L_tot
  float getRatio(const char* specName, std::vector<float> measurements, std::string);

 private:

  float getSpeciesProb(const char* specName, std::vector<float> measurements, std::string gsfClass);
  std::string _name;
  TDirectory * _directory;
  std::vector<LikelihoodSpecies*> _specList;
  std::vector<float> _priorList;
  int _ecalsubdet;
  int _ptbin;

};
#endif
    
