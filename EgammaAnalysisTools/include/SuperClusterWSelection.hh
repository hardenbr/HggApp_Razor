//-----------------------------------------------------
// Description:
// Make a selection of a supercluster and MT to study
// the GSF tracking efficiency with W candle
//-----------------------------------------------------

#ifndef SuperClusterWSelection_HH
#define SuperClusterWSelection_HH

#include "EgammaAnalysisTools/include/EgammaBase.h"
#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "../include/McTruthEvent.hh"

#include <TFile.h>
#include <TTree.h>

class SuperClusterWSelection : public EgammaBase {

public:

  typedef std::pair<unsigned int,unsigned int> aLSSegment;
  typedef std::vector< std::pair<unsigned int,unsigned int> > LSSegments;
  typedef unsigned int aRun;
  typedef std::map< aRun, LSSegments > runsLSSegmentsMap;
  typedef std::pair < aRun, LSSegments > aRunsLSSegmentsMapElement;

  //! constructor
  SuperClusterWSelection(TTree *tree=0);
  //! destructor
  virtual ~SuperClusterWSelection();
  //! set the list of the required triggers
  void requireTrigger(vector<int> requiredTriggers) { m_requiredTriggers = requiredTriggers; }
  //! use a prefix for the output 
  void setPrefix(const char *prefix) { _prefix = prefix; }
  //! set the signal type to select based on MC-truth
  void setSignal(int signal) { m_signal = signal; }
  //! display the efficiency table
  void displayEfficiencies();

  /// Fill RunLSMap according to json file
  void fillRunLSMap();
  /// Set Good Run LS
  void setJsonGoodRunList(const string& jsonFilePath);
  /// check if Run/LS is a good one
  bool isGoodRunLS();

  //! loop over events
  void Loop();

private:

  void ConfigSelection(Selection* selection, Counters* counters);

  //! dump the final tree
  void createOutTree();
  float myMt, myEta, myPt;
  int myRecoGSF;
  TTree *outTree_;
  TFile *fileOut_;
  
  //! the prefix to add to the output files
  const char *_prefix;

  //! the list of required triggers
  vector<int> m_requiredTriggers;

  //! global selection and counters
  Selection *_selection;
  Counters *_counters;

  ///goodRUN/LS list
  runsLSSegmentsMap goodRunLS; 
  std::string jsonFile;

  //! flag Data/MC
  bool isData_;

  //! kind of process to be selected with MC truth in case of MC
  enum { wjets=0, zjets=1, wother=2, zother=3, all=4 };
  int m_signal;

  McTruthEvent mcevent;
  int WToENuDecay;

};


#endif
