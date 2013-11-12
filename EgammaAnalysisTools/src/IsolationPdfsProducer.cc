#include <iostream>
#include <string> 
#include <math.h>

#include "TLorentzVector.h"

#include "CommonTools/include/Utils.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "EgammaAnalysisTools/include/IsolationPdfsProducer.hh"

IsolationPdfsProducer::IsolationPdfsProducer(TTree *tree)
  : EgammaBase(tree) {
  
  std::string fileCuts("/afs/cern.ch/user/e/emanuele/scratch0/Likelihood21X/OfflineAnalysis/EgammaAnalysisTools/config/IsolationPdfsProducer/cuts.txt");
  std::string fileSwitches("/afs/cern.ch/user/e/emanuele/scratch0/Likelihood21X/OfflineAnalysis/EgammaAnalysisTools/config/IsolationPdfsProducer/switches.txt");

  m_selection = new Selection(fileCuts,fileSwitches);
  m_selection->addSwitch("requireTrigger");
  m_selection->addSwitch("applyIDOnProbe");
  m_selection->addCut("meeWindow");
  m_selection->addCut("etaEleAcc");
  m_selection->addCut("ptEleAcc");
  m_selection->addCut("relSumPtTracks");
  m_selection->summary();

}

IsolationPdfsProducer::~IsolationPdfsProducer() { }

void IsolationPdfsProducer::LoopZTagAndProbe(const char *treefilesuffix) {

  if(fChain == 0) return;

  bookHistos();

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // trigger
    Utils anaUtils;
    bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

    //    if(!passedHLT) continue;

    // best tag-probe pair = mee closest to Z mass
    float minpull = 100000.;
    float mass = 1000;
    for(int iele1=0; iele1<nEle; iele1++) {
      TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);
      for(int iele2=iele1+1; iele2<nEle; iele2++) {
        TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);

        if( m_selection->getSwitch("etaEleAcc") && 
            (!m_selection->passCut("etaEleAcc",etaEle[iele1]) ||
             !m_selection->passCut("etaEleAcc",etaEle[iele2]) ) ) continue;

        if( m_selection->getSwitch("ptEleAcc") && 
            (!m_selection->passCut("ptEleAcc",electron1.Pt()) ||
             !m_selection->passCut("ptEleAcc",electron2.Pt()) ) ) continue;

        mass = (electron1+electron2).M();
        m_Zmass->Fill(mass);
        float pull=fabs(mass-91.1876);

        if(pull < minpull) {
          minpull = pull;
          electrons[0] = iele1;
          electrons[1] = iele2;
        }

      }

    }

    // start the tag & probe
    if( m_selection->passCut("meeWindow",minpull) ) {
      
      for(int iele=0; iele<2; ++iele) {
        int tag=electrons[0], probe=electrons[1];
        if(iele=1) {
          tag=electrons[1]; 
          probe=electrons[0];
        }

        TLorentzVector probeP4(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);
        TLorentzVector tagP4(pxEle[tag],pyEle[tag],pzEle[tag],energyEle[tag]);
        
        float eta = etaEle[probe];

        /// define ECAL sub-detector region
        int iecal = (fabs( etaEle[probe])<1.479) ? 0 : 1;
        
        // apply the electron ID loose on the tag electron        
        //        float tagIdentified = eleIdCutsEle[probe];
        float tagIdentified = false; // attention! variable no more in the ntuple

        // Tracker isolation
        // on the tag electron...
        
        bool tagIsolated = ( m_selection->passCut("relSumPtTracks",dr04TkSumPtEle[tag]) );

        // ID on the probe electron...
        bool probeIdentified = true;
        if ( m_selection->getSwitch("applyIDOnProbe") ) {
          //          probeIdentified = eleIdCutsEle[probe];
          probeIdentified = false; // attention! variable no more in the ntuple
        }
        
        float TkSumPtRel = dr04TkSumPtEle[probe];
        float EcalSumEtRel = dr04EcalRecHitSumEtEle[probe];
        float HcalSumEtRel = dr04HcalTowerSumEtEle[probe];
        float EcalRecHitsSumEtRel = 0;

        /// fill the electron ID pdfs only if:
        /// the tag is loose isolated and identified
        /// the probe is loose identified
        if( tagIsolated && tagIdentified && probeIdentified ) {

          tkSumPtRel   [iecal] -> Fill ( TkSumPtRel );
          ecalSumEtRel [iecal] -> Fill ( EcalSumEtRel );
          hcalSumEtRel [iecal] -> Fill ( HcalSumEtRel );
          ecalRecHitsSumEtRel [iecal] -> Fill ( EcalRecHitsSumEtRel );

        } // fill histograms

      } // loop over the 2 Z electrons

    } // end tag and probe

  } // loop over events

}

void IsolationPdfsProducer::LoopQCD() {

  if(fChain == 0) return;

  bookHistos();

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // trigger
//     Utils anaUtils;
//     bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

//     if(!passedHLT) continue;

    // fill the PDFs for QCD with all the (identified) reco'ed electrons
    for(int iele=0;iele<nEle;iele++) {

      if( m_selection->getSwitch("etaEleAcc") && 
          ! m_selection->passCut("etaEleAcc",etaEle[iele]) ) continue;

      // attention, variable no more in the ntuple
//       if ( m_selection->getSwitch("applyIDOnProbe") &&
//            ! eleIdCutsEle[iele] ) continue;
      
      TLorentzVector eleP4(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);

      /// define the bins in which can be splitted the PDFs
      int iecal = (fabs( etaEle[iele])<1.479) ? 0 : 1;

      float TkSumPtRel = dr04TkSumPtEle[iele];
      float EcalSumEtRel = dr04EcalRecHitSumEtEle[iele];
      float HcalSumEtRel = dr04HcalTowerSumEtEle[iele];
      float EcalRecHitsSumEtRel = 0;
      
      tkSumPtRel   [iecal] -> Fill ( TkSumPtRel );
      ecalSumEtRel [iecal] -> Fill ( EcalSumEtRel );
      hcalSumEtRel [iecal] -> Fill ( HcalSumEtRel );
      ecalRecHitsSumEtRel [iecal] -> Fill ( EcalRecHitsSumEtRel );

    } // loop on electrons

  } // loop on events
  
}

// make jet PDFs from W+jets 
void IsolationPdfsProducer::LoopWjets() {

  if(fChain == 0) return;

  bookHistos();

  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    // trigger
//     Utils anaUtils;
//     bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);

//     if(!passedHLT) continue;

    bool tauPresence=false;
    
    int mceleindex = -1;
    for(int imc=0; imc<50; imc++) {
      // not only ele from W->enu: there is enu emission (V_ud?) in madgraph
      if ( (fabs(idMc[imc])==11) ) {
        mceleindex=imc;
        break;
      }
      // since the stable particle list is truncated, if there is a tau 
      // not possible to say what happens...
      if ( (fabs(idMc[imc])==15) ) {
        tauPresence=true;
        break;
      }
    }

    if(tauPresence) continue;

    TVector3 mcEle(0,0,0);
    if(mceleindex>-1) mcEle = TVector3(pMc[mceleindex]*cos(phiMc[mceleindex])*sin(thetaMc[mceleindex]), 
                                       pMc[mceleindex]*sin(phiMc[mceleindex])*sin(thetaMc[mceleindex]),
                                       pMc[mceleindex]*cos(thetaMc[mceleindex]));

    // fill the PDFs for jets excluding the real electron in W+jets
    for(int iele=0;iele<nEle;iele++) {

      TLorentzVector eleP4(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);
      
      float deltaR = 1000;
      if(mceleindex>-1) deltaR = eleP4.Vect().DeltaR(mcEle);
      
      if(deltaR<0.3) continue;
      
      if( m_selection->getSwitch("etaEleAcc") && 
          ! m_selection->passCut("etaEleAcc",etaEle[iele]) ) continue;

      // attntion...
//       if ( m_selection->getSwitch("applyIDOnProbe") &&
//            ! eleIdCutsEle[iele] ) continue;
      
      
      /// define the bins in which can be splitted the PDFs
      int iecal = (fabs( etaEle[iele])<1.479) ? 0 : 1;

      float TkSumPtRel = dr04TkSumPtEle[iele];
      float EcalSumEtRel = dr04EcalRecHitSumEtEle[iele];
      float HcalSumEtRel = dr04HcalTowerSumEtEle[iele];
      float EcalRecHitsSumEtRel = 0;
      
      tkSumPtRel   [iecal] -> Fill ( TkSumPtRel );
      ecalSumEtRel [iecal] -> Fill ( EcalSumEtRel );
      hcalSumEtRel [iecal] -> Fill ( HcalSumEtRel );
      ecalRecHitsSumEtRel [iecal] -> Fill ( EcalRecHitsSumEtRel );

    } // loop on electrons

  } // loop on events
  
}

void IsolationPdfsProducer::bookHistos() {

  m_Zmass = new TH1F("Zmass", "Zmass", 260, 0., 130.);
  
  int nbins = 100;
  
  float tkSumPtRelMin  = 0.0;
  float tkSumPtRelMax = 0.5;
  float ecalSumEtRelMin  = 0.0;
  float ecalSumEtRelMax = 0.5;
  float hcalSumEtRelMin  = 0.0;
  float hcalSumEtRelMax = 0.5;
  float ecalRecHitsSumEtRelMin  = 0.0;
  float ecalRecHitsSumEtRelMax = 0.5;

  // booking histos isolation
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  for (int iecal=0; iecal<2; iecal++) {
    
    char histo[200];
    
    sprintf(histo,"tkSumPtRel_%d",iecal);
    tkSumPtRel[iecal]    = new TH1F(histo, histo, nbins, tkSumPtRelMin, tkSumPtRelMax);
    sprintf(histo,"ecalSumEtRel_%d",iecal);
    ecalSumEtRel[iecal]    = new TH1F(histo, histo, nbins, ecalSumEtRelMin, ecalSumEtRelMax);
    sprintf(histo,"hcalSumEtRel_%d",iecal);
    hcalSumEtRel[iecal]    = new TH1F(histo, histo, nbins, hcalSumEtRelMin, hcalSumEtRelMax);
    sprintf(histo,"ecalRecHitsSumEtRel_%d",iecal);
    ecalRecHitsSumEtRel[iecal]    = new TH1F(histo, histo, nbins, ecalRecHitsSumEtRelMin, ecalRecHitsSumEtRelMax);

  }

}

void IsolationPdfsProducer::saveHistos(const char *filename) {

  char outfilename[200];
  sprintf(outfilename,"%sHistograms.root",filename);
  TFile *file = TFile::Open(outfilename,"recreate");

  m_Zmass->Write();

  for (int iecal=0; iecal<2; iecal++) {

    tkSumPtRel[iecal]->Write();
    ecalSumEtRel[iecal]->Write();
    hcalSumEtRel[iecal]->Write();
    ecalRecHitsSumEtRel[iecal]->Write();

  }

  file->Close();
  
}
