#ifndef ElectronLikelihood_H
#define ElectronLikelihood_H

#include "EgammaAnalysisTools/include/LikelihoodSwitches.h"
#include "EgammaAnalysisTools/include/LikelihoodMeasurements.h"
#include "EgammaAnalysisTools/include/LikelihoodPdfProduct.h"
#include <TDirectory.h>
#include <vector>


class ElectronLikelihood {

 public:
  
  //! ctor, not used for this algo (need initialization from ES)
  ElectronLikelihood () {} ;

  //! ctor
  ElectronLikelihood (TDirectory *EB0lt15dir, TDirectory *EB1lt15dir, TDirectory *EElt15dir,
                      TDirectory *EB0gt15dir, TDirectory *EB1gt15dir, TDirectory *EEgt15dir,
		      LikelihoodSwitches eleIDSwitches,
		      std::string signalWeightSplitting,
		      std::string backgroundWeightSplitting,
		      bool splitSignalPdfs,
		      bool splitBackgroundPdfs) ;

  //! dtor
  virtual ~ElectronLikelihood () ;

  //! get the result of the algorithm
  float result (const LikelihoodMeasurements electron) const;
  //! get the log-expanded result of the algorithm
  float resultLog (const LikelihoodMeasurements electron) const;

 private:

  //! build the likelihood model from histograms 
  //! in Barrel file and Endcap file
  void Setup (TDirectory *EB0lt15dir, TDirectory *EB1lt15dir, TDirectory *EElt15dir,
              TDirectory *EB0gt15dir, TDirectory *EB1gt15dir, TDirectory *EEgt15dir,
	      std::string signalWeightSplitting,
	      std::string backgroundWeightSplitting,
	      bool splitSignalPdfs,
	      bool splitBackgroundPdfs) ;


  //! get the input variables from the electron and the e-Setup
  void getInputVar (std::vector<float> &measuremnts) const ;

  //! likelihood below 15GeV/c
  LikelihoodPdfProduct *_EB0lt15lh, *_EB1lt15lh, *_EElt15lh;
  //! likelihood above 15GeV/c
  LikelihoodPdfProduct *_EB0gt15lh, *_EB1gt15lh, *_EEgt15lh;

  //! general parameters of all the ele id algorithms
  LikelihoodSwitches m_eleIDSwitches ;

  //! splitting rule for PDF's
  std::string m_signalWeightSplitting;
  std::string m_backgroundWeightSplitting;
  bool m_splitSignalPdfs;
  bool m_splitBackgroundPdfs;

};

#endif // ElectronLikelihood_H
