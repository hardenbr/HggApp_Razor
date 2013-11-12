///////////////////////////////////////////////////////////////////////
//
//    FILE: PUWeight.h
//   CLASS: PUWeight
// AUTHORS: I. Gonzalez Caballero
//    DATE: 09/03/2011
//
///////////////////////////////////////////////////////////////////////
#include "CommonTools/include/PUWeight.h"

// ROOT Includes
#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"

// C++ includes
#include <iostream>
using namespace std;

//Set DEBUGPUWEIGHT to 1 to get some debug information. Set it to 2 for more
//detail debug information.
#define DEBUGPUWEIGHT 2

#ifdef DEBUG
#define DEBUGPUWEIGHT 1
#endif

PUWeight::PUWeight(float luminosity, const char* year):
  fData(0),
  fMC(0),
  fWeight(0) {
  
  std::cout << "PU reweighting (old constructor): reading ideal MC" << std::endl;
  
  // Load Data histogram
  if (!LoadDataHistogram(luminosity, year)) return;
  
  // No MC given. Take ideal MC
  IdealMCHistogram();
  
  //Calculate Weight
  CalculateWeight();
}

PUWeight::PUWeight(const char* mcfolder, const char* mcproccess, float luminosity, const char* year):
		   
  fData(0),
  fMC(0),
  fWeight(0) {

  std::cout << "PU reweighting (new constructor): reading both data and MC from file" << std::endl;

  // Load Data histogram
  if (!LoadDataHistogram(luminosity, year)) return;
  
  // Load MC Histogram
  if (!LoadMCHistogram(mcfolder, mcproccess)) return;
  
  // Calculate Weight
  CalculateWeight();
}

TH1F* PUWeight::LoadMCHistogram(const char* mcfolder, const char* mcproccess) {
#ifdef DEBUGPUWEIGHT
  cout << ">> Getting pileup for the MC " << mcproccess << " inside " << mcfolder << "..." << endl;
#endif
  
  TString dsfile;
  //  dsfile.Form("/afs/cern.ch/user/e/emanuele/public/puReweighting/%s/%s_PU.root", mcfolder, mcproccess); // for Summer 11 MC sample
  dsfile.Form("/afs/cern.ch/user/m/mwlebour/public/s6MCPileUp.root"); // for Fall 11 MC sample
      
#if (DEBUGPUWEIGHT > 1)
  cout << "   + Opening " << dsfile << endl;
#endif
  
  TFile* fds = TFile::Open(dsfile);
  if (!fds) {
    cerr << "ERROR [PUWeight]: Could not open file " << dsfile << "!"  << endl
	 << "                  Revise dataset name (" << mcproccess 
	 << ") or internet connection" << endl;
    return 0;
  }
  
  // Read dataset histogram...
#if (DEBUGPUWEIGHT > 1)
  cout << "   + Looking for histogram..." << endl;
#endif
  
  fMC = (TH1F*) fds->Get("pileup")->Clone("PU_MC");
  if (!fMC) {
    cerr << "ERROR [PUWeight]: Could not find histogram for dataset " << mcproccess << "!"
	 << endl;
    return 0;
  }
  fMC->SetDirectory(0);

  if (fMC->Integral() != 1) {
    cout << "NOTE [PUWeight]: MC histogram is not normalized to 1! Normalizing..."
	 << endl;
    fMC->Scale(1./fMC->Integral());
  }

  fds->Close();
  return fMC;  
}

void PUWeight::SetMCHistogram(const TH1F* mcHisto) {

  fMC = (TH1F*)mcHisto->Clone();
  
  if (fMC->Integral() != 1) {
    cout << "NOTE [PUWeight]: MC histogram is not normalized to 1! Normalizing..."
	 << endl;
    fMC->Scale(1./fMC->Integral());
  }

  CalculateWeight();
}

TH1F* PUWeight::LoadDataHistogram(float luminosity, const char* year) {

#ifdef DEBUGPUWEIGHT
  cout << ">> Getting pileup for the " << luminosity << " pb-1 of data..." << endl;       
#endif
  
  TString dtfile;
  TFile* fdt = 0;
  if (luminosity > 0) {
    //    dtfile.Form("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Pileup_2011_to_173692_LPLumiScale_68mb.root");
    // dtfile.Form("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_175832-177515_PromptReco_JSON.pileupTruth_v2.root");
    dtfile.Form("/afs/cern.ch/user/s/sixie/public/Pileup/targets/PUTarget.Full2011.160404-180252.root");
	  
#if (DEBUGPUWEIGHT > 1)
    cout << "   + Opening " << dtfile << endl;
#endif
    
    fdt = TFile::Open(dtfile);
    if (!fdt) {
      cerr << "NOTE [PUWeight]: Could not find file " << dtfile << "!"  << endl;
      cerr << "                 Trying default PU profile for data..." << endl;
    }
  }

  if (!fdt) {
    //    dtfile="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Pileup_2011_to_173692_LPLumiScale_68mb.root";
    // dtfile="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/PileUp/Cert_175832-177515_PromptReco_JSON.pileupTruth_v2.root";
    dtfile="/afs/cern.ch/user/s/sixie/public/Pileup/targets/PUTarget.Full2011.160404-180252.root";

#if (DEBUGPUWEIGHT > 1)
    cout << "   + Opening " << dtfile << endl;
#endif
    
    fdt = TFile::Open(dtfile);
    if (!fdt) {
      cerr << "ERROR [PUWeight]: Could not find default profile in \"" 
	   << dtfile << "\"!"  << endl
	   << "                  Is your internet connection working?" << endl;
      return 0;
    }
  }
  
  //Read data histogram...
  fData = (TH1F*) fdt->Get("pileup")->Clone("PU_Data");
  if (!fData) {
    cerr << "ERROR [PUWeight]: Could not find histogram for data!" << endl;
    return 0;
  }
  
  fData->SetDirectory(0);
  
  if (fData->Integral() != 1) {
    cout << "NOTE [PUWeight]: Data histogram is not normalized to 1! Normalizing..."
	 << endl;
    fData->Scale(1./fData->Integral());
    fData->GetXaxis()->Set(51,-0.5,50.5);
  }

  fdt->Close();

  return fData;
}


TH1F* PUWeight::CalculateWeight() {
  if (fData && fMC) {
    unsigned int nbins = fData->GetXaxis()->GetNbins();
    float xmin  = fData->GetXaxis()->GetXmin();
    float xmax  = fData->GetXaxis()->GetXmax();
    fWeight = new TH1F("PUWeight", "PU Weight", nbins, xmin, xmax);
    fWeight->SetDirectory(0);
    fWeight->Divide(fData, fMC);
  }
  else {
    cerr << "ERROR [PUWeight]: Something weird happened when trying to calculate the weights."
	 << endl 
	 << "                  I could not find the data and/or mc histograms!"
	 << endl;
  }

  return fWeight;
}

TH1F* PUWeight::IdealMCHistogram() {
  unsigned int nbins = 25;
  float xmin = -0.5;
  float xmax = 24.5;

  if (fData) {
    nbins = fData->GetXaxis()->GetNbins();
    xmin  = fData->GetXaxis()->GetXmin();
    xmax  = fData->GetXaxis()->GetXmax();
  }


  fMC = new TH1F("PU_MC", "PU^{MC} Weight", nbins, xmin, xmax);
  float idealpu[] = {0.0698146584, 0.0698146584, 0.0698146584, 
		     0.0698146584, 0.0698146584, 0.0698146584,
		     0.0698146584, 0.0698146584, 0.0698146584,
		     0.0698146584, 0.0698146584, 0.0630151648,
		     0.0526654164, 0.0402754482, 0.0292988928,
		     0.0194384503, 0.0122016783, 0.007207042,
		     0.004003637,  0.0020278322, 0.0010739954,
		     0.0004595759, 0.0002229748, 0.0001028162,
		     4.58337152809607E-05};
  for (unsigned int i = 0; i < nbins; i++) {
    if (i < 25)
      fMC->Fill(i, idealpu[i]);
    else
      fMC->Fill(i, 0);
  }
  return fMC;
}
