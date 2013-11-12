#ifndef JPsiBasicStudiesEle_h
#define JPsiBasicStudiesEle_h

#include <vector>
#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Monitor.hh"
#include "EgammaAnalysisTools/include/EgammaBase.h"
#include <TVector3.h>
#include <TLorentzVector.h>

class PFElectronSeedingEfficiency : public EgammaBase{
public:
  
  //! constructor
  PFElectronSeedingEfficiency(TTree *tree=0);
  //! destructor
  virtual ~PFElectronSeedingEfficiency();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies(std::string filename);
  //! histos
  void bookHistos();
  void writeHistos();
  
private:

  //! for MC
  bool findMcTree();
  //! for PF studies
  int getBin(float eta, float pt);
  int getBin(float pt);
  //! getPt util
  float getPt(float px, float py) { return sqrt(px*px + py*py); }
  //! counters for the efficiencies display
  Counters* counter;
  //! be verbose during runtime
  bool _verbose;
  //! name of rootfile with dataset
  std::string _datasetName;

  //! MC studies
  int _theGenEle;
  TVector3 trueEle;
  
  //! PF studies
  float thr[150], thrPS[20];

  //! histos
  TH1F *H_genEta, *H_genPhi, *H_genPt; 
  TH1F *H_deltaRmin_recoGenTrack;

  TH1F *H_genEtaForEff_befQuality,  *H_genPtForEff_befQuality; 

  TH1F *H_genEtaForEff_all,  *H_genPtForEff_all; 
  TH1F *H_recoEtaForEff_all, *H_recoPtForEff_all; 

  TH1F *H_step1EoP[9],  *H_step1Pt[9]; 
  TH1F *H_step1Chi2[9], *H_step1Hits[9];
  TH1F *H_step1EoP_etaIncl[3],  *H_step1Pt_etaIncl[3]; 
  TH1F *H_step1Chi2_etaIncl[3];

  TH1F *H_genEtaForEff_step1,  *H_genPtForEff_step1; 
  TH1F *H_recoEtaForEff_step1, *H_recoPtForEff_step1; 

  TH1F *H_genEtaForEff_notStep1,  *H_genPtForEff_notStep1; 
  TH1F *H_recoEtaForEff_notStep1, *H_recoPtForEff_notStep1; 

  TH1F *H_step2Chi2[9],         *H_step2Hits[9],         *H_step2Pt[9]; 
  TH1F *H_step2Chi2_etaIncl[3], *H_step2Hits_etaIncl[9], *H_step2Pt_etaIncl[3];  

  TH1F *H_step2_mva;

  TH1F *H_genEtaForEff_step2,  *H_genPtForEff_step2; 
  TH1F *H_recoEtaForEff_step2, *H_recoPtForEff_step2; 

  TH1F *H_genEtaForEff_step2mva,  *H_genPtForEff_step2mva; 
  TH1F *H_recoEtaForEff_step2mva, *H_recoPtForEff_step2mva; 

  TH1F *H_genEtaForEff_full,  *H_genPtForEff_full; 
  TH1F *H_recoEtaForEff_full, *H_recoPtForEff_full; 
};
#endif
