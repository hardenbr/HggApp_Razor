//---------------------------------------
// Description:
// Create the Likelihood Pdfs from Z->ee Tag and probe
//---------------------------------------

#ifndef LHPDFSPRODUCER_H
#define LHPDFSPRODUCER_H

#include <vector>

#include "TH1F.h"
#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "EgammaAnalysisTools/include/EgammaBase.h"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"

#include <TVector3.h>
#include <TMath.h>
#include <TF2.h>
#include <TLorentzVector.h>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>

class LHPdfsProducer : public EgammaBase {

public:

  typedef std::pair<unsigned int,unsigned int> aLSSegment;
  typedef std::vector< std::pair<unsigned int,unsigned int> > LSSegments;
  typedef unsigned int aRun;
  typedef std::map< aRun, LSSegments > runsLSSegmentsMap;
  typedef std::pair < aRun, LSSegments > aRunsLSSegmentsMapElement;

  //! constructor
  LHPdfsProducer(TTree *tree=0);
  //! destructor
  virtual ~LHPdfsProducer();
  //! loop over events doing signal pdfs with Zee - tag and probe 
  void LoopZTagAndProbe(const char *treefilesuffix);
  //! loop over events doing signal pdfs with Zee - tag and probe. Then match with MC truth from Z 
  void LoopZTagAndProbeForMcTruth(const char *treefilesuffix);
  //! loop over events doing signal pdfs with Zee - MC
  void LoopZ(const char *treefilesuffix);
  //! loop over events doing signal pdfs with Zee - MC
  void LoopZwithMass(const char *treefilesuffix);
  //! loop over events doing bkg pdfs on QCD
  void LoopQCD();
  //! loop over events doing bkg pdfs on QCD using the tag and probe method
  void LoopQCDTagAndProbe(const char *treefilesuffix);
  //! loop over events doing bkg pdfs on W+jets
  void LoopWjets();
  //! loop over events doing bkg pdfs on Z+jets
  void LoopZjets(const char *filename);
  //! save the pdfs in a ROOT file
  void saveHistos(const char *filename);
  //! returns the output of the custom cut electron ID                                                                                
  void isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput);

  //! display the efficiency table 
  void displayEfficiencies(const char *filename);

  /// Fill RunLSMap according to json file 
  void fillRunLSMap();
  /// Set Good Run LS 
  void setJsonGoodRunList(const string& jsonFilePath);
  /// check if Run/LS is a good one 
  bool isGoodRunLS();
  /// reload TriggerMask if necessary (data file is changed). Should be called for each event inside the event loop
  bool reloadTriggerMask(bool newVersion=false);
  /// set the list of required trigger to produce the bitmask
  void setRequiredTriggers(const std::vector<std::string>& reqTriggers);
  //check if the event passed HLT. To be called per event
  bool hasPassedHLT();
  //get the value of the requested bits 
  vector<int> getHLTOutput();

private:
  
  //! book the histograms of the PDFs
  void bookHistos();

  //! configurations
  void configSelection(Selection* selection, Counters* counters);

  //! for json
  int isData_;
  runsLSSegmentsMap goodRunLS;
  std::string jsonFile;
  std::string lastFile;

  //! HLT
  std::vector<std::string> requiredTriggers;
  std::vector<int> m_requiredTriggers;

  //! contains the class-dependent electron ID cuts
  std::vector<Selection*> m_EgammaCutIDSelection;

  //! to evaluate custom eleID 
  CutBasedEleIDSelector EgammaCutBasedID;

  //! the tag and the probe
  int electrons[2];
  int muons[2];

  //! the configurable selection
  Selection *m_selection;
  Counters *m_counters;

  // ---------- monitoring histograms ------------
  TH1F *m_Zmass;

  /// Electrons: not splitted
  /// histo[ecalsubdet][ptbin]
  TH1F *dPhiCaloUnsplitEle[2][2];
  TH1F *dPhiVtxUnsplitEle[2][2];
  TH1F *dEtaUnsplitEle[2][2];
  TH1F *EoPoutUnsplitEle[2][2];
  TH1F *HoEUnsplitEle[2][2];  
  TH1F *shapeFisherUnsplitEle[2][2];  
  TH1F *sigmaIEtaIEtaUnsplitEle[2][2];
  TH1F *sigmaIEtaIPhiUnsplitEle[2][2];
  TH1F *sigmaIPhiIPhiUnsplitEle[2][2];
  TH1F *s1s9UnsplitEle[2][2];
  TH1F *s9s25UnsplitEle[2][2];
//   TH1F *dxyUnsplitEle[2][2];
//   TH1F *dxySigUnsplitEle[2][2];
  
  /// Electrons class-splitted
  /// histo[ecalsubdet][ptbin][class]
  TH1F *dPhiCaloClassEle[2][2][2];
  TH1F *dPhiVtxClassEle[2][2][2];
  TH1F *dEtaClassEle[2][2][2];
  TH1F *EoPoutClassEle[2][2][2];
  TH1F *HoEClassEle[2][2][2];
  TH1F *shapeFisherClassEle[2][2][2];
  TH1F *sigmaIEtaIEtaClassEle[2][2][2];
  TH1F *sigmaIEtaIPhiClassEle[2][2][2];
  TH1F *sigmaIPhiIPhiClassEle[2][2][2];
  TH1F *s1s9ClassEle[2][2][2];
  TH1F *s9s25ClassEle[2][2][2];
//   TH1F *dxyClassEle[2][2][2];
//   TH1F *dxySigClassEle[2][2][2];

  /// Electrons fullclass-splitted
  /// histo[ecalsubdet][ptbin][fullclass]
  TH1F *dPhiCaloFullclassEle[2][2][4];
  TH1F *dPhiVtxFullclassEle[2][2][4];
  TH1F *dEtaFullclassEle[2][2][4];
  TH1F *EoPoutFullclassEle[2][2][4];
  TH1F *HoEFullclassEle[2][2][4];
  TH1F *shapeFisherFullclassEle[2][2][4];
  TH1F *sigmaIEtaIEtaFullclassEle[2][2][4];
  TH1F *sigmaIEtaIPhiFullclassEle[2][2][4];
  TH1F *sigmaIPhiIPhiFullclassEle[2][2][4];
  TH1F *s1s9FullclassEle[2][2][4];
  TH1F *s9s25FullclassEle[2][2][4];
  //  TH1F *dxyFullclassEle[2][2][4];
  //  TH1F *dxySigFullclassEle[2][2][4];

};

#endif
