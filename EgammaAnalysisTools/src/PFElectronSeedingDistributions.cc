#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Utils.hh"
#include "EgammaAnalysisTools/include/PFElectronSeedingDistributions.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/EfficiencyEvaluator.hh"

#include <iostream>
#include <fstream>
#include <string>

#include <TTree.h>

using namespace bits;

PFElectronSeedingDistributions::PFElectronSeedingDistributions(TTree *tree) 
  : EgammaBase(tree) { }

PFElectronSeedingDistributions::~PFElectronSeedingDistributions(){ }

bool PFElectronSeedingDistributions::findMcTree() {
  
  _theGenEle=-1;
  int indE=999;
  for(int imc=0;imc<nMc;imc++) {
    if( indE<999 ) continue;
    float thePt = pMc[imc]*fabs(sin(thetaMc[imc]));
    if ( thePt<2 || thePt>50 ) continue;
    if( abs(idMc[imc])==11 ) {indE = imc; }
  }
  
  if( indE<999 ) {
    _theGenEle = indE;
    trueEle.SetPtThetaPhi(pMc[_theGenEle]*fabs(sin(thetaMc[_theGenEle])),thetaMc[_theGenEle],phiMc[_theGenEle]);
    
  }
  
  return ( indE<999 );
}

void PFElectronSeedingDistributions::Loop() {
  
  _verbose=false;
  if(fChain == 0) return;
  
  // to modify here
  bool checkES = false;
  float dRRecoGenTrackCut = 0.15;
  
  // histos
  bookHistos();
  
  // utils to unpack reco flags
  Utils anaUtils;
  
  // contatori eventi
  int events = 0;
  int montecarlo = 0;
  
  // studio inefficienza tracking
  int montecarloEle         = 0;
  int matchedWithTrack      = 0;
  int genTrackPtGt2         = 0;
  int highQualityTrack      = 0; 
  int matchedWithTrackPreId = 0;

  // contatori per PF seeding, step1
  int passed1a = 0;
  int passed1b = 0;
  int passed1c = 0;
  int changed1 = 0;
  int ourEcalMatch        = 0;
  int florianEcalMatch    = 0;
  int passedStep1beforeES = 0;
  int passedStep1withES   = 0;
  int notPassedStep1      = 0;

  // contatori per PF seeding, step1/2
  int goodRangeStep2 = 0;  

  // contatori per PF seeding, step2
  int passed2a = 0;
  int passed2b = 0;
  int ourMatchStep2 = 0;
  int florianMatchStep2 = 0;
  int step2WithMva  = 0;     

  // final cross-check
  int ourFullPreId = 0;
  int florianFullPreId = 0;

  // read thresholds for PF seeding studies
  std::ifstream ifs("config/thresholdPF/Threshold.dat");
  std::ifstream ifsPS("config/thresholdPF/PSThreshold.dat");
  for (int iy=0;iy<72;iy++) ifs   >> thr[iy];
  for (int iy=0;iy<12;iy++) ifsPS >> thrPS[iy];




  // analysis
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000==0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    events++;
    
    // A) TO STUDY TRACKING EFFICIENCY: we look for tracks
      
    // 1) looking for general tracks 
    for(int theT=0; theT<nTrack; theT++) {

      TVector3 pTrack(pxTrack[theT],pyTrack[theT],pzTrack[theT]);

      if(pTrack.Pt()<2.0) continue;
      
      // SEEDING STUDIES
      
      genTrackPtGt2++;
	
      // looking for preId variables related to the matched track
      int preIdIndex = -1;
      for (int iPreId=0; iPreId<nPFpreId; iPreId++) {
        if (trackIndexPFpreId[iPreId]==theT) { preIdIndex = iPreId; }
      }
      
      if (preIdIndex>-1) { 
        matchedWithTrackPreId++;
	    
        // kf track kinematics
        int trk = trackIndexPFpreId[preIdIndex];
        float trackPt  = pTrack.Pt();
        float trackEta = pTrack.Eta();
	
        // for histos
        int binAverEta = -1;
        if (trackPt<6  && trackPt>0)  { binAverEta = 0; } 
        if (trackPt>=6 && trackPt<12) { binAverEta = 1; } 
        if (trackPt>=12)              { binAverEta = 2; } 
	    
        // to study efficiency: before any cut except match with reco-track	and good quality track
        H_recoEtaForEff_all -> Fill(trackEta);
        H_recoPtForEff_all  -> Fill(trackPt);
	
	
        // C) TO STUDY PF SEEDING: 1st STEP
        int ipteta      = getBin(trackEta,trackPt);
        int ibin        = ipteta*8;
        float chi2cut   = thr[ibin+0];
        float ep_cutmin = thr[ibin+1];
        
        // cout << "eta = " << trackEta << ", pt = " << trackPt 
        // << ", chi2chu = " << chi2cut << ", epcut = " << ep_cutmin << endl;  	    
        
        // basic selection for step1
        bool passing1a = chi2MatchPFpreId[preIdIndex]<chi2cut;
        bool passing1b = eopMatchPFpreId[preIdIndex]>ep_cutmin;
        bool passing1c = kfNHitsPFpreId[preIdIndex]>10;
        if(passing1a) passed1a++;
        if(passing1b) passed1b++;
        if(passing1c) passed1c++;
	
        // cross-check
        bool goodMatchStep1Before = passing1a && passing1b && passing1c;
        if (goodMatchStep1Before) ourEcalMatch++;
        if (ecalMatchingPFpreId[preIdIndex]) florianEcalMatch++;
	
        // further selection based on high track pt
        bool goodMatchStep1noES = goodMatchStep1Before;
        bool passing1d = trackPt>50;
        bool passing1e = trackPt<2;
        if (passing1d) goodMatchStep1noES = true;
        if (passing1e) goodMatchStep1noES = false;
        bool changedAfter1 = false;
        if (goodMatchStep1noES!=goodMatchStep1Before) changedAfter1 = true;
        if (changedAfter1) changed1++;
	
        // full step 1, before ES selection
        if (goodMatchStep1noES) passedStep1beforeES++;
	    
        // final step 1, with ES selection
        bool goodMatchStep1 = goodMatchStep1noES;
        bool passingES = psMatchingPFpreId[preIdIndex]; 
        if (checkES && !passingES) goodMatchStep1 = false;
        if (goodMatchStep1) passedStep1withES++;
        
        // histos - input distributions
        // bin by bin
        H_step1Hits[ipteta] -> Fill(kfNHitsPFpreId[preIdIndex]);
        H_step1Pt[ipteta]   -> Fill(trackPt);
        if (kfNHitsPFpreId[preIdIndex]>10) {
          H_step1Chi2[ipteta] -> Fill(chi2MatchPFpreId[preIdIndex]);
          H_step1EoP[ipteta]  -> Fill(eopMatchPFpreId[preIdIndex]);
        }
        // averaged in eta
        H_step1Pt_etaIncl[binAverEta] -> Fill(trackPt);
        if (kfNHitsPFpreId[preIdIndex]>10) {
          H_step1Chi2_etaIncl[binAverEta] -> Fill(chi2MatchPFpreId[preIdIndex]); 
          H_step1EoP_etaIncl[binAverEta]  -> Fill(eopMatchPFpreId[preIdIndex]);
        }
	
        // to study efficiency: after step1
        if (goodMatchStep1) {
          H_recoEtaForEff_step1 -> Fill(trackEta);
          H_recoPtForEff_step1  -> Fill(trackPt);
        }
        if (!goodMatchStep1) {
          H_recoEtaForEff_notStep1 -> Fill(trackEta);
          H_recoPtForEff_notStep1  -> Fill(trackPt);
        }
      
        // D) TO STUDY PF SEEDING TRACKER BASED EFFICIENCY: 2nd step
        bool GoodRange = ((fabs(trackEta)<2.4) && (trackPt>2.));	       
	
        // basic for second step
        int hit1max         = int(thr[ibin+2]);
        float chiredmin     = thr[ibin+3];
        bool passing2a      = kfChi2PFpreId[preIdIndex]>chiredmin;
        bool passing2b      = kfNHitsPFpreId[preIdIndex]<hit1max;	      
        bool goodMatchStep2 = (passing2a || passing2b);
	
        bool passedMva = false;
        if(!goodMatchStep1) { 
          notPassedStep1++;
	  
          if(GoodRange) { 
            goodRangeStep2++;
            
            // real step 2
            if (passing2a) passed2a++;
            if (passing2b) passed2b++;
            if (goodMatchStep2) ourMatchStep2++;
            // cross-check
            if (trackFilteredPFpreId[preIdIndex]) florianMatchStep2++;	
            
            // mva cut
            if (goodMatchStep2) {
              float BDTcut=thr[ibin+4];
              // if (mvaPFpreId[preIdIndex]>BDTcut) passedMva = true;
              if (preidedPFpreId[preIdIndex]) passedMva = true;
            }
            if (passedMva) step2WithMva++;
            
            // histos - input distributions
            // bin by bin
            H_step2Hits[ipteta] -> Fill(kfNHitsPFpreId[preIdIndex]);
            H_step2Pt[ipteta]   -> Fill(trackPt);
            H_step2Chi2[ipteta] -> Fill(kfChi2PFpreId[preIdIndex]);
            // averaged in eta
            H_step2Hits_etaIncl[binAverEta] -> Fill(kfNHitsPFpreId[preIdIndex]);
            H_step2Pt_etaIncl[binAverEta]   -> Fill(trackPt);
            H_step2Chi2_etaIncl[binAverEta] -> Fill(kfChi2PFpreId[preIdIndex]); 
            // mva
            if (goodMatchStep2) {
              //H_step2_mva -> Fill(mvaPFpreId[preIdIndex]);
              H_step2_mva -> Fill(-999);
            }
            
            // to study efficiency: after step2
            if (goodMatchStep2) {
              H_recoEtaForEff_step2 -> Fill(trackEta);
              H_recoPtForEff_step2  -> Fill(trackPt);
              
              // to study the EE+ EE- efficiency asymmetry of the MVA cut
              int iEE=-1;
              if(trackEta<-1.479) iEE=0;
              else if(trackEta>1.479) iEE=1;
              
              if(iEE>-1) {
                H_Chi2EEBeforeMVA[iEE] -> Fill(chi2MatchPFpreId[preIdIndex]);
                H_DeltaEtaEEBeforeMVA[iEE] -> Fill(deltaEtaMatchPFpreId[preIdIndex]);
                H_DeltaPhiEEBeforeMVA[iEE] -> Fill(deltaPhiMatchPFpreId[preIdIndex]);
                H_EoPEEBeforeMVA[iEE] -> Fill(eopMatchPFpreId[preIdIndex]);
              }

              if (passedMva) {
                H_recoEtaForEff_step2mva -> Fill(trackEta);
                H_recoPtForEff_step2mva  -> Fill(trackPt);
              }
            }
          }
        }
	
        // summary
        bool ourFull     = (passedMva || goodMatchStep1);
        bool florianFull = preidedPFpreId[preIdIndex];
        if (ourFull)     ourFullPreId++;
        if (florianFull) florianFullPreId++;
	
        if (ourFull){ 
          H_recoEtaForEff_full -> Fill(trackEta);
          H_recoPtForEff_full  -> Fill(trackPt);
        }
      } // preId found in the tree
    } // track loop
  } // loop over entries
  
  cout << endl;
  cout << endl;
  cout << "eventi totali = " << events << endl; 
  cout << "eventi totali after MC thrut = " << montecarlo << endl; 
  cout << endl;
  cout << "per studiare tracking based PF seeding" << endl;
  cout << "elettroni totali generati (after MC thrut cut) = "   << montecarloEle         << endl;
  cout << "elettroni gene with match with tracker = "           << matchedWithTrack      << endl;
  cout << "elettroni totali generati with pt>2 = "              << genTrackPtGt2         << endl; 
  cout << "la traccia e' high quality = "                       << highQualityTrack      << endl;
  cout << "elettroni gene with match with tracker and preID = " << matchedWithTrackPreId << endl;
  cout << endl;
  cout << "step1 seeding"           << endl;
  cout << "chi2 cut = "             << passed1a << endl;
  cout << "E/P cut = "              << passed1b << endl;
  cout << "nhits cut = "            << passed1c << endl;
  cout << "changed for track pt = " << changed1 << endl;
  cout << endl;
  cout << "our ECAL match before change and ES = "                    << ourEcalMatch        << endl;
  cout << "cross-check: Florian's ECAL match before change and ES = " << florianEcalMatch    << endl;
  cout << "step1 beforeES, after changes = "                          << passedStep1beforeES << endl; 
  cout << "final step1 (after ECAL matching, changes and ES) = "      << passedStep1withES   << endl;         
  cout << endl;
  cout << "not passed step 1 = " << notPassedStep1 << endl;   
  cout << "good range tracks = " << goodRangeStep2 << endl;
  cout << "chi2 cut = "          << passed2a       << endl;
  cout << "nhits cut = "         << passed2b       << endl;
  cout << "our passed step2 = "  << ourMatchStep2  << endl;
  cout << "cross-check: Florian's passed step2 = " << florianMatchStep2 << endl;
  cout << "passed full 2nd step with pre-id = "    << step2WithMva << endl;
  cout << endl;
  cout << "our summary: passed pre-id = " << ourFullPreId << endl;
  cout << "cross-check: Florian's summary: passed pre-id = " << florianFullPreId << endl;
  cout << endl;
  cout << endl;

  writeHistos();
}
  
void PFElectronSeedingDistributions::displayEfficiencies(std::string datasetName) {
  
  std::string::size_type loc = datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    datasetName.erase(loc);
  }
}

void PFElectronSeedingDistributions::bookHistos(){

  char title[500];

  float xbins[17];
  for(int i=0; i<=10; i++) xbins[i]=2+i;
  for(int i=0; i<4; i++) xbins[11+i]=14+2*i;
  xbins[15] = 25; xbins[16] = 30;

  // efficiencies: before any cut except track quality
  H_recoEtaForEff_all = new TH1F("H_recoEtaForEff_all", "eta",  50,-2.5, 2.5);
  H_recoPtForEff_all  = new TH1F("H_recoPtForEff_all",  "p_{T}",16, xbins);
  
  // step1 variables
  for (int ipteta=0; ipteta<9; ipteta++){  
    sprintf(title,"H_step1Chi2[%d]",ipteta);
    H_step1Chi2[ipteta] = new TH1F(title, title, 50, 0.0, 25.); 
    sprintf(title,"H_step1EoP[%d]",ipteta);    
    H_step1EoP[ipteta]  = new TH1F(title, title, 50, 0.3, 1.8);
    sprintf(title,"H_step1Hits[%d]",ipteta);        
    H_step1Hits[ipteta] = new TH1F(title, title, 25, 0.0, 25.);
    sprintf(title,"H_step1Pt[%d]",ipteta);            
    H_step1Pt[ipteta]   = new TH1F(title, title, 50, 0.0, 25.);
  }
  for (int ipteta=0; ipteta<3; ipteta++){  
    sprintf(title,"H_step1Chi2_etaIncl[%d]",ipteta);
    H_step1Chi2_etaIncl[ipteta] = new TH1F(title, title, 50, 0.0, 25.); 
    sprintf(title,"H_step1EoP_etaIncl[%d]",ipteta);    
    H_step1EoP_etaIncl[ipteta]  = new TH1F(title, title, 50, 0.3, 1.8);
    sprintf(title,"H_step1Pt_etaIncl[%d]",ipteta);            
    H_step1Pt_etaIncl[ipteta]   = new TH1F(title, title, 50, 0.0, 25.);
  }

  // efficiencies: after step 1
  H_recoEtaForEff_step1 = new TH1F("H_recoEtaForEff_step1", "eta",  50,-2.5, 2.5);
  H_recoPtForEff_step1  = new TH1F("H_recoPtForEff_step1",  "p_{T}",16, xbins);

  H_recoEtaForEff_notStep1 = new TH1F("H_recoEtaForEff_notStep1", "eta",  50,-2.5, 2.5);
  H_recoPtForEff_notStep1  = new TH1F("H_recoPtForEff_notStep1",  "p_{T}", 16, xbins);

  // step2 variables
  for (int ipteta=0; ipteta<9; ipteta++){  
    sprintf(title,"H_step2Hits[%d]",ipteta);
    H_step2Hits[ipteta] = new TH1F(title, title, 20, 3., 23.);
    sprintf(title,"H_step2Chi2[%d]",ipteta);
    H_step2Chi2[ipteta] = new TH1F(title, title, 50, 0., 10.);
    sprintf(title,"H_step2Pt[%d]",ipteta);            
    H_step2Pt[ipteta]   = new TH1F(title, title, 50, 0.0, 25.);
  }
  for (int ipteta=0; ipteta<3; ipteta++){  
    sprintf(title,"H_step2Chi2_etaIncl[%d]",ipteta);
    H_step2Chi2_etaIncl[ipteta] = new TH1F(title, title, 50, 0.0, 10.); 
    sprintf(title,"H_step2Hits_etaIncl[%d]",ipteta);
    H_step2Hits_etaIncl[ipteta] = new TH1F(title, title, 20, 3., 23.);
    sprintf(title,"H_step2Pt_etaIncl[%d]",ipteta);            
    H_step2Pt_etaIncl[ipteta]   = new TH1F(title, title, 50, 0.0, 25.);
  }
  // id
  H_step2_mva = new TH1F("H_step2_mva", "H_step2_mva", 50, -1., +1.);

  // to check the EE asymmetry before MVA cut
  for (int ieta=0; ieta<2; ieta++){  
    sprintf(title,"H_Chi2EEBeforeMVA[%d]",ieta);
    H_Chi2EEBeforeMVA[ieta] = new TH1F(title, title, 100, 0.0, 50.); 
    sprintf(title,"H_DeltaEtaEEBeforeMVA[%d]",ieta);
    H_DeltaEtaEEBeforeMVA[ieta] = new TH1F(title, title, 100, 0.0, 2*TMath::Pi()); 
    sprintf(title,"H_DeltaPhiEEBeforeMVA[%d]",ieta);
    H_DeltaPhiEEBeforeMVA[ieta] = new TH1F(title, title, 100, 0.0, 2*TMath::Pi()); 
    sprintf(title,"H_EoPEEBeforeMVA[%d]",ieta);    
    H_EoPEEBeforeMVA[ieta]  = new TH1F(title, title, 100, 0.3, 3.0);
  }

  
  // efficiencies: after step 2 without mva
  H_recoEtaForEff_step2 = new TH1F("H_recoEtaForEff_step2", "eta",  50,-2.5, 2.5);
  H_recoPtForEff_step2  = new TH1F("H_recoPtForEff_step2",  "p_{T}", 16, xbins);

  // efficiencies: after step 2 with mva
  H_recoEtaForEff_step2mva = new TH1F("H_recoEtaForEff_step2mva", "eta",  50,-2.5, 2.5);
  H_recoPtForEff_step2mva  = new TH1F("H_recoPtForEff_step2mva",  "p_{T}", 16, xbins);

  // efficiencies: full pre-id
  H_recoEtaForEff_full = new TH1F("H_recoEtaForEff_full", "eta",  50,-2.5, 2.5);
  H_recoPtForEff_full  = new TH1F("H_recoPtForEff_full",  "p_{T}", 16, xbins);
}

void PFElectronSeedingDistributions::writeHistos(){

  TFile *fOut = new TFile("histos.root","RECREATE");
  
  for (int ipteta=0; ipteta<9; ipteta++){  
    H_step1Chi2[ipteta] -> Write();
    H_step1EoP[ipteta]  -> Write();
    H_step1Hits[ipteta] -> Write();
    H_step1Pt[ipteta]   -> Write();
  }
  
  for (int ipteta=0; ipteta<3; ipteta++){  
    H_step1Chi2_etaIncl[ipteta] -> Write();
    H_step1EoP_etaIncl[ipteta]  -> Write();
    H_step1Pt_etaIncl[ipteta]   -> Write();
  }
  
  for (int ipteta=0; ipteta<9; ipteta++){  
    H_step2Chi2[ipteta] -> Write();
    H_step2Hits[ipteta] -> Write();
    H_step2Pt[ipteta]   -> Write();
  }
  
  for (int ipteta=0; ipteta<3; ipteta++){  
    H_step2Chi2_etaIncl[ipteta] -> Write();
    H_step2Hits_etaIncl[ipteta] -> Write();
    H_step2Pt_etaIncl[ipteta]   -> Write();
  }
  
  H_step2_mva -> Write();
  
  for(int ieta=0; ieta<2; ieta++) {
    H_Chi2EEBeforeMVA[ieta]     -> Write();
    H_DeltaEtaEEBeforeMVA[ieta] -> Write();
    H_DeltaPhiEEBeforeMVA[ieta] -> Write();
    H_EoPEEBeforeMVA[ieta]      -> Write();
  }

  // efficiency evaluation
  EfficiencyEvaluator effEta1("efficiencyEta1.root","RECREATE");
  effEta1.AddNumerator(H_recoEtaForEff_all);
  effEta1.AddNumerator(H_recoEtaForEff_step1);
  effEta1.SetDenominator(H_recoEtaForEff_all);
  effEta1.ComputeEfficiencies();
  effEta1.SetTitle("PF preID efficiency vs #eta - step1");
  effEta1.SetXaxisTitle("track #eta");
  effEta1.SetYaxisTitle("efficiency");
  effEta1.SetYaxisMin(0.0);
  effEta1.Write();
  
  EfficiencyEvaluator effEta2("efficiencyEta2.root","RECREATE");
  effEta2.AddNumerator(H_recoEtaForEff_all);
  effEta2.AddNumerator(H_recoEtaForEff_notStep1);
  effEta2.AddNumerator(H_recoEtaForEff_step2);
  effEta2.AddNumerator(H_recoEtaForEff_step2mva);
  effEta2.SetDenominator(H_recoEtaForEff_all);
  effEta2.ComputeEfficiencies();
  effEta2.SetTitle("PF preID efficiency vs #eta - step2");
  effEta2.SetXaxisTitle("track #eta");
  effEta2.SetYaxisTitle("efficiency");
  effEta2.SetYaxisMin(0.0);
  effEta2.Write();
  
  EfficiencyEvaluator effEta3("efficiencyEta3.root","RECREATE");
  effEta3.AddNumerator(H_recoEtaForEff_all);
  effEta3.AddNumerator(H_recoEtaForEff_full);
  effEta3.SetDenominator(H_recoEtaForEff_all);
  effEta3.ComputeEfficiencies();
  effEta3.SetTitle("PF preID efficiency vs #eta - step or step2");
  effEta3.SetXaxisTitle("track #eta");
  effEta3.SetYaxisTitle("efficiency");
  effEta3.SetYaxisMin(0.0);
  effEta3.Write();
  
  EfficiencyEvaluator effPt1("efficiencyPt1.root","RECREATE");
  effPt1.AddNumerator(H_recoPtForEff_all);
  effPt1.AddNumerator(H_recoPtForEff_step1);
  effPt1.SetDenominator(H_recoPtForEff_all);
  effPt1.ComputeEfficiencies();
  effPt1.SetTitle("PF preID efficiency vs p_{T} - step1");
  effPt1.SetXaxisTitle("track p_{T} [GeV]");
  effPt1.SetYaxisTitle("efficiency");
  effPt1.SetYaxisMin(0.0);
  effPt1.Write();
  
  EfficiencyEvaluator effPt2("efficiencyPt2.root","RECREATE");
  effPt2.AddNumerator(H_recoPtForEff_all);
  effPt2.AddNumerator(H_recoPtForEff_notStep1);
  effPt2.AddNumerator(H_recoPtForEff_step2);
  effPt2.AddNumerator(H_recoPtForEff_step2mva);
  effPt2.SetDenominator(H_recoPtForEff_all);
  effPt2.ComputeEfficiencies();
  effPt2.SetTitle("PF preID efficiency vs p_{T} - step2");
  effPt2.SetXaxisTitle("track p_{T}");
  effPt2.SetYaxisTitle("efficiency");
  effPt2.SetYaxisMin(0.0);
  effPt2.Write();
  
  EfficiencyEvaluator effPt3("efficiencyPt3.root","RECREATE");
  effPt3.AddNumerator(H_recoPtForEff_all);
  effPt3.AddNumerator(H_recoPtForEff_full);
  effPt3.SetDenominator(H_recoPtForEff_all);
  effPt3.ComputeEfficiencies();
  effPt3.SetTitle("PF preID efficiency vs p_{T} - step or step2");
  effPt3.SetXaxisTitle("track p_{T}");
  effPt3.SetYaxisTitle("efficiency");
  effPt3.SetYaxisMin(0.0);
  effPt3.Write();
}

int PFElectronSeedingDistributions::getBin(float eta, float pt){
  int ie=0;
  int ip=0;
  if (fabs(eta)<1.2) ie=0;
  else{ if (fabs(eta)<1.68) ie=1;
  else ie=2;
  }
  if (pt<6) ip=0;
  else {  if (pt<12) ip=1;     
  else ip=2;
  }
  int iep= ie*3+ip;
  return iep;
}

int PFElectronSeedingDistributions::getBin(float pt){
  int ip=0;
  if (pt<6) ip=0;
  else {  if (pt<12) ip=1;
  else ip=2;
  }
  return ip;
}
