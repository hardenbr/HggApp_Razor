#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>

ElectronLikelihood::ElectronLikelihood (TDirectory *EB0lt15dir, TDirectory *EB1lt15dir, TDirectory *EElt15dir,
                                        TDirectory *EB0gt15dir, TDirectory *EB1gt15dir, TDirectory *EEgt15dir,
					LikelihoodSwitches eleIDSwitches,
					std::string signalWeightSplitting,
					std::string backgroundWeightSplitting,
					bool splitSignalPdfs,
					bool splitBackgroundPdfs) :
  _EB0lt15lh (new LikelihoodPdfProduct ("electronID_EB0_ptLt15_likelihood",0,0)) ,
  _EB1lt15lh (new LikelihoodPdfProduct ("electronID_EB1_ptLt15_likelihood",1,0)) ,
  _EElt15lh (new LikelihoodPdfProduct ("electronID_EE_ptLt15_likelihood",2,0)) ,
  _EB0gt15lh (new LikelihoodPdfProduct ("electronID_EB0_ptGt15_likelihood",0,1)) ,
  _EB1gt15lh (new LikelihoodPdfProduct ("electronID_EB1_ptGt15_likelihood",1,1)) ,
  _EEgt15lh (new LikelihoodPdfProduct ("electronID_EE_ptGt15_likelihood",2,1)) ,
  m_eleIDSwitches (eleIDSwitches) ,
  m_signalWeightSplitting (signalWeightSplitting), 
  m_backgroundWeightSplitting (backgroundWeightSplitting),
  m_splitSignalPdfs (splitSignalPdfs), 
  m_splitBackgroundPdfs (splitBackgroundPdfs)  
{
  Setup (EB0lt15dir, EB1lt15dir, EElt15dir,
         EB0gt15dir, EB1gt15dir, EEgt15dir,
	 signalWeightSplitting, backgroundWeightSplitting,
	 splitSignalPdfs, splitBackgroundPdfs) ;
}



// --------------------------------------------------------



ElectronLikelihood::~ElectronLikelihood () {
  delete _EB0lt15lh ;
  delete _EB1lt15lh ;
  delete _EElt15lh ;
  delete _EB0gt15lh ;
  delete _EB1gt15lh ;
  delete _EEgt15lh ;
}



// --------------------------------------------------------


void 
ElectronLikelihood::Setup (TDirectory *EB0lt15dir, TDirectory *EB1lt15dir, TDirectory *EElt15dir,
                           TDirectory *EB0gt15dir, TDirectory *EB1gt15dir, TDirectory *EEgt15dir,
			   std::string signalWeightSplitting,
			   std::string backgroundWeightSplitting,
			   bool splitSignalPdfs,
			   bool splitBackgroundPdfs) 
{

  // ECAL BARREL0 (|eta|<0.8) LIKELIHOOD - Pt < 15 GeV region
  _EB0lt15lh->initFromFile (EB0lt15dir) ;

  _EB0lt15lh->addSpecies ("electrons", 1.0) ;
  _EB0lt15lh->addSpecies ("hadrons",  1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EB0lt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EB0lt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EB0lt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EB0lt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EB0lt15lh->setSplitFrac ("electrons", "class0", 1.0) ;
    _EB0lt15lh->setSplitFrac ("electrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EB0lt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EB0lt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EB0lt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EB0lt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EB0lt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EB0lt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EB0lt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EB0lt15lh->addPdf ("electrons", "OneOverEMinusOneOverP",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EB0lt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EB0lt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EB0lt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EB0lt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EB0lt15lh->setSplitFrac ("hadrons", "class0", 1.0) ;
    _EB0lt15lh->setSplitFrac ("hadrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EB0lt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EB0lt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EB0lt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EB0lt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EB0lt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EB0lt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EB0lt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EB0lt15lh->addPdf ("hadrons", "OneOverEMinusOneOverP",         splitBackgroundPdfs) ;



  // ECAL BARREL0 (|eta|<0.8) LIKELIHOOD - Pt >= 15 GeV region
  _EB0gt15lh->initFromFile (EB0gt15dir) ;

  _EB0gt15lh->addSpecies ("electrons", 1.0 ) ;  
  _EB0gt15lh->addSpecies ("hadrons",1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EB0gt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EB0gt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EB0gt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EB0gt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EB0gt15lh->setSplitFrac ("electrons", "class0", 1.0);
    _EB0gt15lh->setSplitFrac ("electrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EB0gt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EB0gt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EB0gt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EB0gt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EB0gt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EB0gt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EB0gt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EB0gt15lh->addPdf ("electrons", "OneOverEMinusOneOverP",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EB0gt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EB0gt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EB0gt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EB0gt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EB0gt15lh->setSplitFrac ("hadrons", "class0", 1.0);
    _EB0gt15lh->setSplitFrac ("hadrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EB0gt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EB0gt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EB0gt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EB0gt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EB0gt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EB0gt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EB0gt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EB0gt15lh->addPdf ("hadrons", "OneOverEMinusOneOverP",         splitBackgroundPdfs) ;


  // ECAL BARREL1 (|eta|>0.8) LIKELIHOOD - Pt < 15 GeV region
  _EB1lt15lh->initFromFile (EB1lt15dir) ;

  _EB1lt15lh->addSpecies ("electrons", 1.0) ;
  _EB1lt15lh->addSpecies ("hadrons",  1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EB1lt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EB1lt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EB1lt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EB1lt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EB1lt15lh->setSplitFrac ("electrons", "class0", 1.0) ;
    _EB1lt15lh->setSplitFrac ("electrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EB1lt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EB1lt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EB1lt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EB1lt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EB1lt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EB1lt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EB1lt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EB1lt15lh->addPdf ("electrons", "OneOverEMinusOneOverP",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EB1lt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EB1lt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EB1lt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EB1lt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EB1lt15lh->setSplitFrac ("hadrons", "class0", 1.0) ;
    _EB1lt15lh->setSplitFrac ("hadrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EB1lt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EB1lt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EB1lt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EB1lt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EB1lt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EB1lt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EB1lt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EB1lt15lh->addPdf ("hadrons", "OneOverEMinusOneOverP",         splitBackgroundPdfs) ;



  // ECAL BARREL1 (|eta|>0.8) LIKELIHOOD - Pt >= 15 GeV region
  _EB1gt15lh->initFromFile (EB1gt15dir) ;

  _EB1gt15lh->addSpecies ("electrons", 1.0 ) ;  
  _EB1gt15lh->addSpecies ("hadrons",1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EB1gt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EB1gt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EB1gt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EB1gt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EB1gt15lh->setSplitFrac ("electrons", "class0", 1.0);
    _EB1gt15lh->setSplitFrac ("electrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EB1gt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EB1gt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EB1gt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EB1gt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EB1gt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EB1gt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EB1gt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EB1gt15lh->addPdf ("electrons", "OneOverEMinusOneOverP",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EB1gt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EB1gt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EB1gt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EB1gt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EB1gt15lh->setSplitFrac ("hadrons", "class0", 1.0);
    _EB1gt15lh->setSplitFrac ("hadrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EB1gt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EB1gt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EB1gt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EB1gt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EB1gt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EB1gt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EB1gt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EB1gt15lh->addPdf ("hadrons", "OneOverEMinusOneOverP",         splitBackgroundPdfs) ;


  // ECAL ENDCAP LIKELIHOOD - Pt < 15 GeV
  _EElt15lh->initFromFile (EElt15dir) ;

  _EElt15lh->addSpecies ("electrons", 1.0 ) ;
  _EElt15lh->addSpecies ("hadrons",1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EElt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EElt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EElt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EElt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EElt15lh->setSplitFrac ("electrons", "class0", 1.0);
    _EElt15lh->setSplitFrac ("electrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now";
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EElt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EElt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EElt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EElt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EElt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EElt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EElt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EElt15lh->addPdf ("electrons", "OneOverEMinusOneOverP",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EElt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EElt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EElt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EElt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EElt15lh->setSplitFrac ("hadrons", "class0", 1.0);
    _EElt15lh->setSplitFrac ("hadrons", "class1", 1.0);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EElt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EElt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EElt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EElt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EElt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EElt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EElt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EElt15lh->addPdf ("hadrons", "OneOverEMinusOneOverP",         splitBackgroundPdfs) ;

  // ECAL ENDCAP LIKELIHOOD - Pt >= 15 GeV
  _EEgt15lh->initFromFile (EEgt15dir) ;

  _EEgt15lh->addSpecies ("electrons", 1.0 ) ;
  _EEgt15lh->addSpecies ("hadrons",1.0) ;

  if(signalWeightSplitting.compare("fullclass")==0) {
    _EEgt15lh->setSplitFrac ("electrons", "fullclass0", 1.0) ;
    _EEgt15lh->setSplitFrac ("electrons", "fullclass1", 1.0) ;
    _EEgt15lh->setSplitFrac ("electrons", "fullclass2", 1.0) ;
    _EEgt15lh->setSplitFrac ("electrons", "fullclass3", 1.0) ;
  }
  else if(signalWeightSplitting.compare("class")==0) {
    _EEgt15lh->setSplitFrac ("electrons", "class0", 1.0) ;
    _EEgt15lh->setSplitFrac ("electrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EEgt15lh->addPdf ("electrons", "dPhi",          splitSignalPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EEgt15lh->addPdf ("electrons", "dEta",          splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EEgt15lh->addPdf ("electrons", "EoP",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EEgt15lh->addPdf ("electrons", "HoE",           splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EEgt15lh->addPdf ("electrons", "sigmaIEtaIEta", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EEgt15lh->addPdf ("electrons", "sigmaIPhiIPhi", splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EEgt15lh->addPdf ("electrons", "fBrem",         splitSignalPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EEgt15lh->addPdf ("electrons", "OneOverEMinusOneOverP",         splitSignalPdfs) ;

  if(backgroundWeightSplitting.compare("fullclass")==0) {
    _EEgt15lh->setSplitFrac ("hadrons", "fullclass0", 1.0) ;
    _EEgt15lh->setSplitFrac ("hadrons", "fullclass1", 1.0) ;
    _EEgt15lh->setSplitFrac ("hadrons", "fullclass2", 1.0) ;
    _EEgt15lh->setSplitFrac ("hadrons", "fullclass3", 1.0) ;
  }
  else if(backgroundWeightSplitting.compare("class")==0) {
    _EEgt15lh->setSplitFrac ("hadrons", "class0", 1.0) ;
    _EEgt15lh->setSplitFrac ("hadrons", "class1", 1.0) ;
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  if (m_eleIDSwitches.m_useDeltaPhi)     _EEgt15lh->addPdf ("hadrons", "dPhi",          splitBackgroundPdfs) ;   
  if (m_eleIDSwitches.m_useDeltaEta)     _EEgt15lh->addPdf ("hadrons", "dEta",          splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useEoverP)       _EEgt15lh->addPdf ("hadrons", "EoP",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useHoverE)       _EEgt15lh->addPdf ("hadrons", "HoE",           splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaEtaEta)  _EEgt15lh->addPdf ("hadrons", "sigmaIEtaIEta", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useSigmaPhiPhi)  _EEgt15lh->addPdf ("hadrons", "sigmaIPhiIPhi", splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useFBrem)        _EEgt15lh->addPdf ("hadrons", "fBrem",         splitBackgroundPdfs) ;
  if (m_eleIDSwitches.m_useOneOverEMinusOneOverP)        _EEgt15lh->addPdf ("hadrons", "OneOverEMinusOneOverP",         splitBackgroundPdfs) ;

}



// --------------------------------------------------------




// --------------------------------------------------------



float 
ElectronLikelihood::result (const LikelihoodMeasurements electron) const 
{

  //=======================================================
  // used classification:
  // nbrem clusters = 0         =>  0
  // nbrem clusters >= 1        =>  1
  //=======================================================

  std::vector<float> measurements ;
  if(m_eleIDSwitches.m_useDeltaPhi) measurements.push_back( electron.deltaPhi );
  if(m_eleIDSwitches.m_useDeltaEta) measurements.push_back( electron.deltaEta );
  if(m_eleIDSwitches.m_useEoverP) measurements.push_back( electron.eSuperClusterOverP );
  if(m_eleIDSwitches.m_useHoverE) measurements.push_back( electron.hadronicOverEm );
  if(m_eleIDSwitches.m_useSigmaEtaEta) measurements.push_back( electron.sigmaIEtaIEta );
  if(m_eleIDSwitches.m_useSigmaPhiPhi) measurements.push_back( electron.sigmaIPhiIPhi );
  if(m_eleIDSwitches.m_useFBrem) measurements.push_back( electron.fBrem );
  if(m_eleIDSwitches.m_useOneOverEMinusOneOverP) measurements.push_back( electron.OneOverEMinusOneOverP );

  // Split using only the 1 / >1 cluster
  int nBremClusters=electron.nBremClusters;
  int bitVal = (nBremClusters==0) ? 0 : 1 ;
  
  char className[20] ;
  if(m_signalWeightSplitting.compare("class")==0) {
    sprintf (className,"class%d",bitVal);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  int subdet = electron.subdet;
  float thisPt =  electron.pt;

  if (subdet==0 && thisPt<15.)
    return _EB0lt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==0 && thisPt>=15.)
    return _EB0gt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==1 && thisPt<15.)
    return _EB1lt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==1 && thisPt>=15.)
    return _EB1gt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==2 && thisPt<15.)
    return _EElt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==2 && thisPt>=15.)
    return _EEgt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else return -999. ;
}

float 
ElectronLikelihood::resultLog (const LikelihoodMeasurements electron) const 
{

  //=======================================================
  // used classification:
  // nbrem clusters = 0         =>  0
  // nbrem clusters >= 1        =>  1
  //=======================================================

  std::vector<float> measurements ;
  if(m_eleIDSwitches.m_useDeltaPhi) measurements.push_back( electron.deltaPhi );
  if(m_eleIDSwitches.m_useDeltaEta) measurements.push_back( electron.deltaEta );
  if(m_eleIDSwitches.m_useEoverP) measurements.push_back( electron.eSuperClusterOverP );
  if(m_eleIDSwitches.m_useHoverE) measurements.push_back( electron.hadronicOverEm );
  if(m_eleIDSwitches.m_useSigmaEtaEta) measurements.push_back( electron.sigmaIEtaIEta );
  if(m_eleIDSwitches.m_useSigmaPhiPhi) measurements.push_back( electron.sigmaIPhiIPhi );
  if(m_eleIDSwitches.m_useFBrem) measurements.push_back( electron.fBrem );
  if(m_eleIDSwitches.m_useOneOverEMinusOneOverP) measurements.push_back( electron.OneOverEMinusOneOverP );

  // Split using only the 1 / >1 cluster
  int nBremClusters=electron.nBremClusters;
  int bitVal = (nBremClusters==0) ? 0 : 1 ;
  
  char className[20] ;
  if(m_signalWeightSplitting.compare("class")==0) {
    sprintf (className,"class%d",bitVal);
  }
  else {
    std::cout << "Only class (non-showering / showering)"
              << " and fullclass (golden / bigbrem / narrow / showering)" 
              << " splitting is implemented right now" << std::endl;
    exit(0);
  }

  int subdet = electron.subdet;
  float thisPt =  electron.pt;

  float lh=-999.;

  if (subdet==0 && thisPt<15.)
    lh = _EB0lt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==0 && thisPt>=15.)
    lh = _EB0gt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==1 && thisPt<15.)
    lh = _EB1lt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==1 && thisPt>=15.)
    lh = _EB1gt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==2 && thisPt<15.)
    lh = _EElt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else if (subdet==2 && thisPt>=15.)
    lh = _EEgt15lh->getRatio ("electrons",measurements,std::string (className)) ;
  else lh = -999. ;

  if(lh<=0) return -20.;
  else if(lh==1) return 20.;
  else return log(lh/(1.0-lh));
  
}
