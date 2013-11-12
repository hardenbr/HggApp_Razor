//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed H->WW->2e2nu
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
// Original code:
//    CutAnaHiggs_2e2nu.cpp
//-------------------------------------------------------

// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>

// Offline analysis includes
#include "EgammaAnalysisTools/include/Application.hh"
#include "CommonTools/include/TriggerMask.hh"

#if Application == 1
#include "EgammaAnalysisTools/include/LikelihoodAnalysis.hh"
#endif
#if Application == 2
#include "EgammaAnalysisTools/include/LHPdfsProducer.hh"
#endif
#if Application == 3
#include "EgammaAnalysisTools/include/sPlotsPdfsComparison.hh"
#endif
#if Application == 4
#include "EgammaAnalysisTools/include/SuperClusterWSelection.hh"
#endif
#if Application == 5
#include "EgammaAnalysisTools/include/sPlotsPdfsComparison.h"
#endif
#if Application == 6
#include "EgammaAnalysisTools/MLFit/src/prepareQCDhistos.h"
#endif
#if Application == 8
#include "EgammaAnalysisTools/include/ZeeTagAndProbe.hh"
#endif
#if Application == 9
#include "EgammaAnalysisTools/include/TestTurnOnCurve.hh"
#endif
#if Application == 10
#include "EgammaAnalysisTools/include/FakeElectronSelector.hh"
#endif
#if Application == 11
#include "EgammaAnalysisTools/include/FakeElectronSelectorWenuPlusOneJet.hh"
#include "EgammaAnalysisTools/include/FakeElectronSelectorWmunuPlusOneJet.hh"
#endif
#if Application == 12
#include "EgammaAnalysisTools/include/FakeElectronSelectorZllPlusOneFake.hh"
#endif
#if Application == 13
#include "EgammaAnalysisTools/include/HZZ4LElectronSelector.hh"
#endif
#if Application == 14
#include "EgammaAnalysisTools/include/TestTurnOnCurveGamma.hh"
#endif

int main(int argc, char* argv[]) {

  char inputFileName[300];
  char outputFileName[300];

#if Application != 5
  if ( argc < 4 ){
    std::cout << "missing argument: insert at least the inputFile with list of root files"    << std::endl; 
    std::cout << "EgammaAnalysis inputFile outputFile 0" << std::endl;
    return 1;
  }
  strcpy(inputFileName, argv[1]);
  strcpy(outputFileName,argv[2]);
  int signal   = atoi(argv[3]);

  // -------------------------
  // loading file:
  TChain *theChain = 0;

#if Application != 3
  theChain = new TChain("ntp1");
#else 
  theChain = new TChain(outputFileName);
#endif

  char Buffer[500];
  char MyRootFile[2000];  
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);
  // get the tree with the conditions from the first file
  // TTree *treeCond = new TTree();
  // int nfiles=1;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer)))
      {
	sscanf(Buffer,"%s",MyRootFile);
	// theChain->Add("rfio:"+TString(MyRootFile));
	theChain->Add(TString(MyRootFile));
	std::cout << "chaining " << MyRootFile << std::endl;
	// if ( nfiles==1 ) {
	//  TFile *firstfile = TFile::Open(MyRootFile);
	//  treeCond = (TTree*)firstfile->Get("Conditions");
	//}
	//nfiles++;
      }
  }
  inputFile->close();
  delete inputFile;
#endif

#if Application == 1

  LikelihoodAnalysis analysis(theChain);
  // analysis.reproduceEgammaCutID();
  // analysis.findEquivalentLHCut( 0.79935 );        // tight eleID
  // analysis.findEquivalentLHCut( 0.957 );          // loose eleID
  // analysis.estimateIDEfficiency(outputFileName);
  // analysis.estimateFakeRate(outputFileName);
  //analysis.estimateFakeRateQCD(outputFileName); 
  analysis.estimateFakeRateForHToWW_EGAMMA(outputFileName);
#endif

#if Application == 2

  char title[1000];
  
  LHPdfsProducer producer(theChain);
  
  // TriggerMask maskSignal(treeCond);
  // maskSignal.requireTrigger("HLT_Ele15_SW_L1R");
  // maskSignal.requireTrigger("HLT_Ele20_SW_L1R");
  // std::vector<int> requiredSignalTriggers = maskSignal.getBits();
  // producer.requireSignalTrigger(requiredSignalTriggers);

  std::vector<std::string> mask;
  mask.push_back("HLT_Jet15U");   
  // mask.push_back("HLT_Jet30U");   
  // mask.push_back("HLT_Jet50U");   
  producer.setRequiredTriggers(mask);
  
  //  sprintf(title,"%s_zTandP_tree.root",outputFileName);  
  // producer.LoopZTagAndProbe(title);
  // sprintf(title,"%s_zTandP_histos.root",outputFileName);    
  // producer.saveHistos(title);
  
  // sprintf(title,"%s_zMC_tree.root",outputFileName);  
  // producer.LoopZTagAndProbeForMcTruth(title);
  // sprintf(title,"%s_zMC_histos.root",outputFileName);    
  // producer.saveHistos(title);

  sprintf(title,"%s_qcdTandP_tree.root",outputFileName);  
  producer.LoopQCDTagAndProbe(title);
  sprintf(title,"%s_qcdTandP_counters.root",outputFileName);
  producer.displayEfficiencies(title);
  sprintf(title,"%s_qcdTandP_histos.root",outputFileName);    
  producer.saveHistos(title);

#endif

#if Application == 3

  sPlotsPdfsComparison p(theChain);
  p.RunOverMC(false);
  p.Loop();

#endif

#if Application == 4

  SuperClusterWSelection p(theChain);
  p.setPrefix(outputFileName);
  p.setSignal(signal);
  p.Loop();
  p.displayEfficiencies();

#endif

#if Application == 5

  int doSignalsPlots = 1;
  int isMC = 1;
  int typeClass = 2;
  int isTP = 0;
  int randomizeNotUsedPdf=0;
  
  TFile *fileData, *fileMC;

  if(doSignalsPlots) {
    fileData = TFile::Open((std::string("results_data/sPlots/Wenu_tree.root")).c_str());
    fileMC = TFile::Open((std::string("wcandleresults/treesW/WjetsMadgraph.root")).c_str());
  } else {
    fileData = TFile::Open((std::string("results_data/sPlots/Wenu_bkgFit_tree.root")).c_str());
    fileMC = TFile::Open((std::string("results/treesW/QCD_Wenu.root")).c_str());    
  }

  TTree *tree = 0;
  if(isMC) tree = (TTree*) fileMC->Get("T1");
  else tree = (TTree*) fileData->Get("dataset");

  std::cout << "tree " << tree << std::endl;
  if (!tree)
    {
      std::cout << "Tree not found" << std::endl;
      exit(-1);
    }

  sPlotsPdfsComparison p;
  //  p.Init(tree, isMC, isTP, randomizeNotUsedPdf, typeClass);
  p.Init(tree, isMC, isTP);
  p.doSignalsPlots(doSignalsPlots);
  p.Loop();

#endif

#if Application == 6

  int isMC = 1;
  
  TFile *fileData, *fileMC;
  fileData = TFile::Open((std::string("results_data/mergedTree.root")).c_str());
  fileMC = TFile::Open((std::string("results/trees/QCD_Pt-20_TuneD6T_tree.root")).c_str());

  TTree *tree = 0;
  if (isMC) tree = (TTree*) fileMC->Get("T1");
  else tree = (TTree*) fileData->Get("T1");

  prepareQCDhistos p;
  p.Init(tree);
  p.Loop();

#endif

#if Application == 8

  ZeeTagAndProbe analysis(theChain);
  analysis.Loop(outputFileName);

#endif

#if Application == 9

  TestTurnOnCurve analysis(theChain);
  analysis.measureTurnOn(outputFileName);
  
#endif

#if Application == 10

  FakeElectronSelector analysis(theChain);
  analysis.Loop(outputFileName);

#endif

#if Application == 11

  TString basename(outputFileName);

  std::cout << "Selecting W->enu + fake electrons..." << std::endl;
  FakeElectronSelectorWenuPlusOneJet wenuanalysis(theChain);
  TString wenuj = basename + TString("-Wenu-");
  wenuanalysis.Loop(wenuj.Data());

  std::cout << "Selecting W->munu + fake electrons..." << std::endl;
  FakeElectronSelectorWmunuPlusOneJet wmunuanalysis(theChain);
  TString wmunuj = basename + TString("-Wmunu-");
  wmunuanalysis.Loop(wmunuj);

#endif

#if Application == 12

  FakeElectronSelectorZllPlusOneFake analysis(theChain);
  analysis.Loop(outputFileName);

#endif

#if Application == 13

  HZZ4LElectronSelector analysis(theChain);
  analysis.Loop(outputFileName);

#endif

#if Application == 14

  TestTurnOnCurveGamma analysis(theChain);
  analysis.measureTurnOnGamma(outputFileName);
  
#endif

  return 0;

}
