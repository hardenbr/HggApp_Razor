#define sPlotsPdfsComparison_cxx
#include "include/sPlotsPdfsComparison.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TRandom3.h>

#include <iostream>
#include <math.h>

using namespace std;

void sPlotsPdfsComparison::Loop()
{
  //   In a ROOT session, you can do:
  //      Root > .L sPlotsPdfsComparison.C
  //      Root > sPlotsPdfsComparison t
  //      Root > t.GetEntry(12); // Fill t data members with entry number 12
  //      Root > t.Show();       // Show values of entry 12
  //      Root > t.Show(16);     // Read and show values of entry 16
  //      Root > t.Loop();       // Loop on all entries
  //
  
  if (fChain == 0) return;

  bookHistosFixedBinning();
  bookFullHistos();  

  // configuring electron likelihood
  TFile *fileLH = TFile::Open("pdfs_new.root");
  TDirectory *LHdir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useOneOverEMinusOneOverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useHoverE = false;
  LH = new ElectronLikelihood(&(*LHdir), &(*LHdir), &(*LHdir), &(*LHdir), &(*LHdir), &(*LHdir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);

  Long64_t nentries = fChain->GetEntries();
  //  Long64_t nentries = 200000;
  
  int firstEvent = 0;
  int lastEvent  = nentries;
  cout << firstEvent << " " << lastEvent << endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=firstEvent; jentry<lastEvent;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%10000 == 0) std::cout << ">>> Processing event # " << jentry << " / " << nentries << std::endl;

    // da cambiare (solo x dati, selezionare sgn o fondo)
    //    float weight = (m_doSignal) ? N_sig_sw : N_qcd_sw;
    float wgt = (m_isMC) ? f_weight : weight;
    //    float wgt = 1.0;
    //    if(m_isDataZTaP) wgt = 1.0;
    // float wgt = (m_isMC) ? 1.0 : N_qcd_sw;

    int jecal, jptbin, jclass;
    if(!m_isDataZTaP) { // W->enu sample
      if(m_isMC) {
        if(wgt>500) continue;
        if(fabs(f_eta)<1.0) jecal = 0;
        else if (fabs(f_eta)<1.442) jecal = 1;
        else if (fabs(f_eta)>1.56) jecal = 2;
	else jecal = -1;
        if(f_pt>=10 && f_pt<=20) jptbin = 0;
        else if(f_pt>20) jptbin = 1;
        else continue;
        
        //        if(f_pt<25 || f_pt>45) continue;

        if(f_nbrem==0) jclass = 0;
        else jclass = 1;
//         if(m_doSignal) {
//           if(f_nJets>0 || f_pfmet/f_pt<0.3 || f_isIsoWP80<0) continue;
//         }
      } else {
        if(see<1e-4) continue; // remove residual spikes (if the cleaning is not done)
        if(fabs(f_eta)<1.0) jecal = 0;
        else if (fabs(f_eta)<1.442) jecal = 1;
        else if (fabs(f_eta)>1.56) jecal = 2;
        else jecal = -1;
        if(pt>=10 && pt<=20) jptbin = 0;
        else if(pt>20) jptbin = 1;
        else continue;
        if(nbrem==0) jclass = 0;
        else jclass = 1;
        if(m_doSignal) {
          //          if(nJets>0 || pfmet/pt<0.3 || isIsoWP80<0) continue;
        }
      }
      if(m_isDataZTaP) {
        if(fabs(ztap_eta)<1.0) jecal = 0;
        else if (fabs(ztap_eta)<1.442) jecal = 1;
        else if (fabs(ztap_eta)>1.56) jecal = 2;
        else jecal = -1;
      }

      if (jecal == -1)
	continue;
      // populate the unusued PDFs with uniform distribution (including the under/over-flows)
//       float minSPP = sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinLowEdge(0)-1.0E-03;
//       float maxSPP = sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinUpEdge(sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetNbinsX()+1)+1.0E-03;
//       float minFBrem = fBremClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinLowEdge(0)-1.0E-03;
//       float maxFBrem = fBremClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinUpEdge(fBremClassEle [jecal][jptbin][jclass]->GetNbinsX()+1)+1.0E-03;
      //under and over flows are taken into account in bin 0 and nbins+1
      float minSPP = sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinLowEdge(0);
      float maxSPP = sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinUpEdge(sigmaIPhiIPhiClassEle [jecal][jptbin][jclass]->GetNbinsX()+1);
      float minFBrem = fBremClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinLowEdge(0);
      float maxFBrem = fBremClassEle [jecal][jptbin][jclass]->GetXaxis()->GetBinUpEdge(fBremClassEle [jecal][jptbin][jclass]->GetNbinsX()+1);
      TRandom3 rnd(jentry);

      if(m_isMC) {
        dPhiEle          [jecal] -> Fill ( f_dphi, wgt );
        dEtaEle          [jecal] -> Fill ( f_deta, wgt );
        EoPEle           [jecal] -> Fill ( f_eop, wgt );
        OneOverEMinusOneOverPEle[jecal] -> Fill ( 1./(f_eop*f_pin) - 1./f_pin, wgt );
        HoEEle           [jecal] -> Fill ( f_hoe, wgt );
        sigmaIEtaIEtaEle [jecal] -> Fill ( f_see, wgt );
        fbremEle         [jecal] -> Fill ( f_fbrem, wgt );
        etaEle                   -> Fill ( f_eta, wgt );
        phiEle           [jecal] -> Fill ( f_phi, wgt );
        chargeEle        [jecal] -> Fill ( f_charge, wgt );

        // PDFs for the likelihood
        dPhiUnsplitEle [jecal][jptbin] -> Fill ( f_dphi, wgt );
        dEtaUnsplitEle [jecal][jptbin] -> Fill ( f_deta, wgt );
        EoPUnsplitEle [jecal][jptbin] -> Fill ( f_eop, wgt );
	OneOverEMinusOneOverPUnsplitEle [jecal][jptbin] -> Fill ( 1./(f_eop*f_pin) - 1./f_pin , wgt );
	HoEUnsplitEle  [jecal][jptbin] -> Fill ( f_hoe, wgt );
        sigmaIEtaIEtaUnsplitEle [jecal][jptbin] -> Fill ( f_see, wgt );
        sigmaIPhiIPhiUnsplitEle [jecal][jptbin] -> Fill ( f_spp, wgt ); 
        fBremUnsplitEle [jecal][jptbin] -> Fill ( f_fbrem, wgt );


        dPhiClassEle [jecal][jptbin][jclass] -> Fill ( f_dphi, wgt );
        dEtaClassEle [jecal][jptbin][jclass] -> Fill ( f_deta, wgt );
        EoPClassEle [jecal][jptbin][jclass] -> Fill ( f_eop, wgt );
	
        OneOverEMinusOneOverPClassEle [jecal][jptbin][jclass] -> Fill ( 1./(f_eop*f_pin) - 1./f_pin , wgt );

        HoEClassEle  [jecal][jptbin][jclass] -> Fill ( f_hoe, wgt );
        sigmaIEtaIEtaClassEle [jecal][jptbin][jclass] -> Fill ( f_see, wgt );
        
        if(jclass==0) {
          sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( f_spp, wgt ); 
          fBremClassEle [jecal][jptbin][jclass] -> Fill ( rnd.Uniform(minFBrem,maxFBrem) );
        } else {
	  //  	  if (jecal!=0)
	  sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( rnd.Uniform(minSPP,maxSPP) );
	  //  	  else
	  // 	    sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( f_spp, wgt ); 
	  if (jecal!=0)
	    fBremClassEle [jecal][jptbin][jclass] -> Fill ( f_fbrem, wgt );
	  else
	    fBremClassEle [jecal][jptbin][jclass] -> Fill ( rnd.Uniform(minFBrem,maxFBrem) );
        }
      } else {
        dPhiEle          [jecal] -> Fill ( dphi, wgt );
        dEtaEle          [jecal] -> Fill ( deta, wgt );
        EoPEle           [jecal] -> Fill ( eop, wgt );
	OneOverEMinusOneOverPEle[jecal] -> Fill ( 1./(eop*pin) - 1./pin, wgt );
        HoEEle           [jecal] -> Fill ( hoe, wgt );
        sigmaIEtaIEtaEle [jecal] -> Fill ( see, wgt );
        fbremEle         [jecal] -> Fill ( fbrem, wgt );
        etaEle                   -> Fill ( eta, wgt );
        phiEle           [jecal] -> Fill ( phi, wgt );
        chargeEle        [jecal] -> Fill ( charge, wgt );

        // PDFs for the likelihood
        dPhiUnsplitEle [jecal][jptbin] -> Fill ( dphi, wgt );
        dEtaUnsplitEle [jecal][jptbin] -> Fill ( deta, wgt );
        EoPUnsplitEle [jecal][jptbin] -> Fill ( eop, wgt );
	OneOverEMinusOneOverPUnsplitEle [jecal][jptbin] -> Fill ( 1./(eop*pin) - 1./pin , wgt );
        HoEUnsplitEle  [jecal][jptbin] -> Fill ( hoe, wgt );
        sigmaIEtaIEtaUnsplitEle [jecal][jptbin] -> Fill ( see, wgt );
        sigmaIPhiIPhiUnsplitEle [jecal][jptbin] -> Fill ( spp, wgt );
        fBremUnsplitEle [jecal][jptbin] -> Fill ( fbrem, wgt );

        dPhiClassEle [jecal][jptbin][jclass] -> Fill ( dphi, wgt );
        dEtaClassEle [jecal][jptbin][jclass] -> Fill ( deta, wgt );
        EoPClassEle [jecal][jptbin][jclass] -> Fill ( eop, wgt );
        OneOverEMinusOneOverPClassEle [jecal][jptbin][jclass] -> Fill ( 1./(eop*pin) - 1./pin , wgt );
        HoEClassEle  [jecal][jptbin][jclass] -> Fill ( hoe, wgt );
        sigmaIEtaIEtaClassEle [jecal][jptbin][jclass] -> Fill ( see, wgt );
        if(jclass==0) {
          sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( spp, wgt ); 
          fBremClassEle [jecal][jptbin][jclass] -> Fill ( rnd.Uniform(minFBrem,maxFBrem) );
        } else {
          sigmaIPhiIPhiClassEle [jecal][jptbin][jclass] -> Fill ( rnd.Uniform(minSPP,maxSPP) );
          fBremClassEle [jecal][jptbin][jclass] -> Fill ( fbrem, wgt );
        }
      }
      // the likelihood output
      float f_lh = likelihoodRatio(m_isMC,*LH);
      lhEle [jecal] -> Fill ( f_lh, wgt );
      lhUnsplitEle [jecal][jptbin] -> Fill ( f_lh, wgt );
      lhClassEle [jecal][jptbin][jclass] -> Fill ( f_lh, wgt );
    } else {
      if(fabs(ztap_eta)<1.0) jecal = 0;
      else if (fabs(ztap_eta)<1.442) jecal = 1;
      else jecal = 2;
      dPhiEle          [jecal] -> Fill ( ztap_dphi, wgt );
      dEtaEle          [jecal] -> Fill ( ztap_deta, wgt );
      EoPEle           [jecal] -> Fill ( ztap_eop, wgt );
      OneOverEMinusOneOverPEle[jecal] -> Fill ( 1./(ztap_eop*ztap_pin) - 1./ztap_pin, wgt );
      HoEEle           [jecal] -> Fill ( ztap_hoe, wgt );
      sigmaIEtaIEtaEle [jecal] -> Fill ( ztap_see, wgt );
      fbremEle         [jecal] -> Fill ( ztap_fbrem, wgt );
      etaEle                   -> Fill ( ztap_eta, wgt );
      phiEle           [jecal] -> Fill ( 0.0, wgt );
      chargeEle        [jecal] -> Fill ( ztap_charge, wgt );
    }
  } // loop over events

  char buffer[200];
  char bufferLikelihood[200];
  char hypothesis[200];
  if(m_doSignal) sprintf(hypothesis,"electrons");
  else sprintf(hypothesis,"hadrons");
  if(! m_isDataZTaP ) { // W sample
    if( m_isMC ) {
      sprintf(buffer,"pdfs_%s_histograms_MC.root",hypothesis);
      sprintf(bufferLikelihood,"pdfs_%s_MC.root",hypothesis);
    }
    else {
      sprintf(buffer,"pdfs_%s_histograms_data.root",hypothesis);
      sprintf(bufferLikelihood,"pdfs_%s_data.root",hypothesis);
    }
  } else {
    sprintf(buffer,"pdfs_%s_histograms_ZTaP.root",hypothesis);
  }

  // file for simple plots
  TFile *fileOut = TFile::Open(buffer,"recreate");

  etaEle->Write();
  
  for (int jecal=0; jecal<3; jecal++) {
    dPhiEle[jecal]->Write();
    dEtaEle[jecal]->Write();
    EoPEle[jecal]->Write();           
    OneOverEMinusOneOverPEle[jecal]->Write();           
    HoEEle[jecal]->Write();
    sigmaIEtaIEtaEle[jecal]->Write();
    fbremEle[jecal]->Write();
    phiEle[jecal]->Write();    
    chargeEle[jecal]->Write();    
    lhEle[jecal]->Write();
  }
  
  fileOut->Close();


  // ROOT file with Electron Likelihood PDFs
  TFile *fileOut2 = TFile::Open(bufferLikelihood, "recreate");
  
  for (int iecal=0; iecal<3; iecal++) {

    for(int iptbin=0; iptbin<2; iptbin++) {
      
      dPhiUnsplitEle[iecal][iptbin]->Write();
      dEtaUnsplitEle[iecal][iptbin]->Write();
      EoPUnsplitEle[iecal][iptbin]->Write();
      OneOverEMinusOneOverPUnsplitEle[iecal][iptbin]->Write();
      HoEUnsplitEle[iecal][iptbin]->Write();
      sigmaIEtaIEtaUnsplitEle[iecal][iptbin]->Write();
      sigmaIPhiIPhiUnsplitEle[iecal][iptbin]->Write();
      fBremUnsplitEle[iecal][iptbin]->Write();
      lhUnsplitEle[iecal][iptbin]->Write();

      for(int iclass=0; iclass<2; iclass++) {
      
	dPhiClassEle[iecal][iptbin][iclass]->Write();
	dEtaClassEle[iecal][iptbin][iclass]->Write();
	EoPClassEle[iecal][iptbin][iclass]->Write();
	OneOverEMinusOneOverPClassEle[iecal][iptbin][iclass]->Write();
	HoEClassEle[iecal][iptbin][iclass]->Write();
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass]->Write();
	sigmaIPhiIPhiClassEle[iecal][iptbin][iclass]->Write();
        fBremClassEle[iecal][iptbin][iclass]->Write();
        lhClassEle[iecal][iptbin][iclass]->Write();

      }
    }
  }

  fileOut2->Close();

}

void sPlotsPdfsComparison::bookHistosVariableBinning() {
  
  float dPhiBins[12] = {-0.1, -0.0348 , -0.0272 , -0.0196 , -0.012 , -0.0044 , 0.0032 , 0.0108 , 0.0184 , 
                        0.026 , 0.0336 , 0.1};

  float dEtaBins[12] = {-0.01 , -0.00424  , -0.00272 , -0.00196 , -0.0012 , -0.00044 , 0.00032 , 0.00108 , 
                        0.00184 , 0.0026 , 0.00412 , 0.01};

  float EoPBinsEB[12] = {0, 0.5, 0.75, 0.875, 1, 1.125, 1.25, 1.5, 2, 4.0, 7., 10.};
  float EoPBinsEE[12] = {0, 0.5, 0.75, 0.875, 1, 1.125, 1.25, 1.5, 2, 4.0, 7., 10.};

  float OneOverEMinusOneOverPBinsEB[12] = {-1., -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3};
  float OneOverEMinusOneOverPBinsEE[12] = {-1., -0.8, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3};

  float seeBinsEB[14] = {0., 0.005 , 0.0055 , 0.006 , 0.0065 , 0.007 , 0.0075 , 0.008 , 0.0085 , 0.009 , 0.0095, 0.01,
                         0.015, 0.020 };

  float seeBinsEE[14] = {0., 0.015 , 0.02, 0.021 , 0.022 , 0.023 , 0.024 , 0.025 , 0.026 , 0.027 , 0.028 , 0.029, 
                         0.03, 0.04};

  float HoEBins[6] = {0., 0.005, 0.01, 0.03, 0.06, 0.1 };

  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  char histo[200];
  for (int iecal=0; iecal<3; iecal++) {
      
    sprintf(histo,"dPhiClass_electrons_%d",iecal);
    dPhiEle[iecal] = new TH1F(histo, histo, 10, dPhiBins);
    sprintf(histo,"dEtaClass_electrons_%d",iecal);
    dEtaEle[iecal] = new TH1F(histo, histo, 10, dEtaBins);
    sprintf(histo,"HoEClass_electrons_%d",iecal);
    HoEEle[iecal] = new TH1F(histo, histo, 5, HoEBins);

  }

  
  int iecal = 0;
  sprintf(histo,"EoPClass_electrons_%d",iecal);
  EoPEle[iecal] = new TH1F(histo, histo, 11, EoPBinsEB);
  sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",iecal);
  OneOverEMinusOneOverPEle[iecal] = new TH1F(histo, histo, 11, OneOverEMinusOneOverPBinsEB);
  sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
  sigmaIEtaIEtaEle[iecal] = new TH1F(histo, histo, 13, seeBinsEB);

  iecal = 1;
  sprintf(histo,"EoPClass_electrons_%d",iecal);
  EoPEle[iecal] = new TH1F(histo, histo, 11, EoPBinsEB);
  sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",iecal);
  OneOverEMinusOneOverPEle[iecal] = new TH1F(histo, histo, 11, OneOverEMinusOneOverPBinsEB);
  sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
  sigmaIEtaIEtaEle[iecal] = new TH1F(histo, histo, 13, seeBinsEB);

  iecal = 2;
  sprintf(histo,"EoPClass_electrons_%d",iecal);
  EoPEle[iecal] = new TH1F(histo, histo, 11, EoPBinsEE);
  sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",iecal);
  OneOverEMinusOneOverPEle[iecal] = new TH1F(histo, histo, 11, OneOverEMinusOneOverPBinsEE);
  sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",iecal);
  sigmaIEtaIEtaEle[iecal] = new TH1F(histo, histo, 13, seeBinsEE);

}



void sPlotsPdfsComparison::bookHistosFixedBinning() {
  
  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  char histo[200];
  //  for (int iecal=0; iecal<2; iecal++) {

  sprintf(histo,"etaClass_electrons");
  etaEle = new TH1F(histo, histo, 30, -3.0, 3.0);
  etaEle->Sumw2();

  if(m_doSignal) {
    // EB histograms |eta|<1.0
    sprintf(histo,"dPhiClass_electrons_%d",0);
    dPhiEle[0] = new TH1F(histo, histo, 50, -0.05, 0.05);
    dPhiEle[0]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",0);
    dEtaEle[0] = new TH1F(histo, histo, 50, -0.008, 0.008);
    dEtaEle[0]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",0);
    HoEEle[0] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEEle[0]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",0);
    EoPEle[0] = new TH1F(histo, histo, 50, 0., 3.);
    EoPEle[0]->Sumw2();
    sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",0);
    OneOverEMinusOneOverPEle[0] = new TH1F(histo, histo, 50, -0.2, 0.1);
    OneOverEMinusOneOverPEle[0]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",0);
    sigmaIEtaIEtaEle[0] = new TH1F(histo, histo, 50, 0., 0.03);
    sigmaIEtaIEtaEle[0]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",0);
    fbremEle[0] = new TH1F(histo, histo, 40, -1.0, 1.0);
    fbremEle[0]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",0);
    phiEle[0] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiEle[0]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",0);
    chargeEle[0] = new TH1F(histo, histo, 2, -2., 2.);
    chargeEle[0]->Sumw2();
    sprintf(histo,"lhClass_electrons_%d",0);
    lhEle[0] = new TH1F(histo, histo, 50, 0., 1.);
    lhEle[0]->Sumw2();

    // EB histograms |eta|>1.0
    sprintf(histo,"dPhiClass_electrons_%d",1);
    dPhiEle[1] = new TH1F(histo, histo, 50, -0.05, 0.05);
    dPhiEle[1]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",1);
    dEtaEle[1] = new TH1F(histo, histo, 50, -0.008, 0.008);
    dEtaEle[1]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",1);
    HoEEle[1] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEEle[1]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",1);
    EoPEle[1] = new TH1F(histo, histo, 50, 0., 3.);
    EoPEle[1]->Sumw2();
    sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",1);
    OneOverEMinusOneOverPEle[1] = new TH1F(histo, histo, 50, -0.2, 0.1);
    OneOverEMinusOneOverPEle[1]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",1);
    sigmaIEtaIEtaEle[1] = new TH1F(histo, histo, 50, 0., 0.03);
    sigmaIEtaIEtaEle[1]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",1);
    fbremEle[1] = new TH1F(histo, histo, 40, -1.0, 1.0);
    fbremEle[1]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",1);
    phiEle[1] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiEle[1]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",1);
    chargeEle[1] = new TH1F(histo, histo, 2, -2., 2.);
    chargeEle[1]->Sumw2();
    sprintf(histo,"lhClass_electrons_%d",1);
    lhEle[1] = new TH1F(histo, histo, 50, 0., 1.);
    lhEle[1]->Sumw2();



    // EE histograms
    sprintf(histo,"dPhiClass_electrons_%d",1);
    dPhiEle[2] = new TH1F(histo, histo, 50, -0.1, 0.1);
    dPhiEle[2]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",1);
    dEtaEle[2] = new TH1F(histo, histo, 50, -0.03, 0.03);
    dEtaEle[2]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",1);
    HoEEle[2] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEEle[2]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",1);
    EoPEle[2] = new TH1F(histo, histo, 50, 0., 5.);
    EoPEle[2]->Sumw2();
    sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",1);
    OneOverEMinusOneOverPEle[2] = new TH1F(histo, histo, 50, -0.2, 0.1);
    OneOverEMinusOneOverPEle[2]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",1);
    sigmaIEtaIEtaEle[2] = new TH1F(histo, histo, 50, 0.01, 0.05);
    sigmaIEtaIEtaEle[2]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",1);
    fbremEle[2] = new TH1F(histo, histo, 30, -1.0, 1.0);
    fbremEle[2]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",1);
    phiEle[2] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiEle[2]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",1);
    chargeEle[2] = new TH1F(histo, histo, 2, -2., 2.);
    chargeEle[2]->Sumw2();
    sprintf(histo,"lhClass_electrons_%d",1);
    lhEle[2] = new TH1F(histo, histo, 50, 0., 1.);
    lhEle[2]->Sumw2();
  } else {
    // EB histograms |eta|<1.0
    sprintf(histo,"dPhiClass_electrons_%d",0);
    dPhiEle[0] = new TH1F(histo, histo, 50, -0.15, 0.15);
    dPhiEle[0]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",0);
    dEtaEle[0] = new TH1F(histo, histo, 50, -0.02, 0.02);
    dEtaEle[0]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",0);
    HoEEle[0] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEEle[0]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",0);
    EoPEle[0] = new TH1F(histo, histo, 50, 0., 10.);
    EoPEle[0]->Sumw2();
    sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",0);
    OneOverEMinusOneOverPEle[0] = new TH1F(histo, histo, 50, -0.2, 0.1);
    OneOverEMinusOneOverPEle[0]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",0);
    sigmaIEtaIEtaEle[0] = new TH1F(histo, histo, 50, 0., 0.03);
    sigmaIEtaIEtaEle[0]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",0);
    fbremEle[0] = new TH1F(histo, histo, 50, -1.0, 1.0);
    fbremEle[0]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",0);
    phiEle[0] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiEle[0]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",0);
    chargeEle[0] = new TH1F(histo, histo, 2, -2., 2.);
    chargeEle[0]->Sumw2();
    sprintf(histo,"lhClass_electrons_%d",0);
    lhEle[0] = new TH1F(histo, histo, 50, 0., 1.);
    lhEle[0]->Sumw2();

    // EB histograms |eta|>=1.0
    sprintf(histo,"dPhiClass_electrons_%d",0);
    dPhiEle[1] = new TH1F(histo, histo, 50, -0.15, 0.15);
    dPhiEle[1]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",0);
    dEtaEle[1] = new TH1F(histo, histo, 50, -0.02, 0.02);
    dEtaEle[1]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",0);
    HoEEle[1] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEEle[1]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",0);
    EoPEle[1] = new TH1F(histo, histo, 50, 0., 10.);
    EoPEle[1]->Sumw2();
    sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",0);
    OneOverEMinusOneOverPEle[1] = new TH1F(histo, histo, 50, -0.2, 0.1);
    OneOverEMinusOneOverPEle[1]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",0);
    sigmaIEtaIEtaEle[1] = new TH1F(histo, histo, 50, 0., 0.03);
    sigmaIEtaIEtaEle[1]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",0);
    fbremEle[1] = new TH1F(histo, histo, 50, -1.0, 1.0);
    fbremEle[1]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",0);
    phiEle[1] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiEle[1]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",0);
    chargeEle[1] = new TH1F(histo, histo, 2, -2., 2.);
    chargeEle[1]->Sumw2();
    sprintf(histo,"lhClass_electrons_%d",0);
    lhEle[1] = new TH1F(histo, histo, 50, 0., 1.);
    lhEle[1]->Sumw2();

    // EE histograms
    sprintf(histo,"dPhiClass_electrons_%d",1);
    dPhiEle[2] = new TH1F(histo, histo, 50, -0.3, 0.3);
    dPhiEle[2]->Sumw2();
    sprintf(histo,"dEtaClass_electrons_%d",1);
    dEtaEle[2] = new TH1F(histo, histo, 50, -0.02, 0.02);
    dEtaEle[2]->Sumw2();
    sprintf(histo,"HoEClass_electrons_%d",1);
    HoEEle[2] = new TH1F(histo, histo, 50, 0., 0.15);
    HoEEle[2]->Sumw2();
    sprintf(histo,"EoPClass_electrons_%d",1);
    EoPEle[2] = new TH1F(histo, histo, 50, 0., 10.);
    EoPEle[2]->Sumw2();
    sprintf(histo,"OneOverEMinusOneOverPClass_electrons_%d",1);
    OneOverEMinusOneOverPEle[2] = new TH1F(histo, histo, 50, -0.2, 0.1);
    OneOverEMinusOneOverPEle[2]->Sumw2();
    sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d",1);
    sigmaIEtaIEtaEle[2] = new TH1F(histo, histo, 50, 0.01, 0.08);
    sigmaIEtaIEtaEle[2]->Sumw2();
    sprintf(histo,"fbremClass_electrons_%d",1);
    fbremEle[2] = new TH1F(histo, histo, 50, -1.0, 1.0);
    fbremEle[2]->Sumw2();
    sprintf(histo,"phiClass_electrons_%d",1);
    phiEle[2] = new TH1F(histo, histo, 50, 0., 2*TMath::Pi());
    phiEle[2]->Sumw2();
    sprintf(histo,"chargeClass_electrons_%d",1);
    chargeEle[2] = new TH1F(histo, histo, 2, -2., 2.);
    chargeEle[2]->Sumw2();
    sprintf(histo,"lhClass_electrons_%d",1);
    lhEle[2] = new TH1F(histo, histo, 50, 0., 1.);
    lhEle[2]->Sumw2();
  }
  
}

void sPlotsPdfsComparison::bookFullHistos() {

  int nbins = 200;

  float dPhiMin  = -0.1; // pixelMatchGsfElectron pre-selection: |dPhi| < 0.1
  float dPhiMax  =  0.1;
  float dEtaMin     = -0.015;
  float dEtaMax     =  0.015;
  float EoPMin   =  0.0;
  float EoPMax   =  10.0;
  float OneOverEMinusOneOverPMin   =  -0.1;
  float OneOverEMinusOneOverPMax   =  0.05;
  float HoEMin      =  0.0; // zero-suppression in HCAL
  float HoEMax      =  0.15; // pixelMatchGsfElectron pre-selection: H/E<0.15
  float sigmaIEtaIEtaEBMin = 0.0;
  float sigmaIEtaIEtaEBMax = 0.035;
  float sigmaIEtaIEtaEEMin = 0.01;
  float sigmaIEtaIEtaEEMax = 0.08;
  float sigmaIPhiIPhiEBMin = 0.0;
  float sigmaIPhiIPhiEBMax = 0.03;
  float sigmaIPhiIPhiEEMin = 0.01;
  float sigmaIPhiIPhiEEMax = 0.09;
  float fBremMin = -0.2;
  float fBremMax = 1.0;
  float lhMin = 0.0;
  float lhMax = 1.0;

  char hypothesis[200];
  if(m_doSignal) sprintf(hypothesis,"electrons");
  else sprintf(hypothesis,"hadrons");
  
  // booking histos eleID
  // iecal = 0 --> barrel |eta|< 1.0
  // iecal = 1 --> barrel |eta|>= 1.0
  // iecal = 2 --> endcap
  for (int iecal=0; iecal<3; iecal++) {

    // iptbin = 0: < 15 GeV
    // iptbin = 1: > 15 GeV
    for(int iptbin=0; iptbin<2; iptbin++) {
      
      char histo[200];

      sprintf(histo,"dPhiUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
      dPhiUnsplitEle[iecal][iptbin]     = new TH1F(histo, histo, nbins, dPhiMin, dPhiMax);
      dPhiUnsplitEle[iecal][iptbin]->Sumw2();
      sprintf(histo,"dEtaUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
      dEtaUnsplitEle[iecal][iptbin]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
      dEtaUnsplitEle[iecal][iptbin]->Sumw2();
      sprintf(histo,"EoPUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
      EoPUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, EoPMin, EoPMax);
      EoPUnsplitEle[iecal][iptbin]->Sumw2();
      sprintf(histo,"OneOverEMinusOneOverPUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
      OneOverEMinusOneOverPUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, OneOverEMinusOneOverPMin, OneOverEMinusOneOverPMax);
      OneOverEMinusOneOverPUnsplitEle[iecal][iptbin]->Sumw2();
      sprintf(histo,"HoEUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
      HoEUnsplitEle[iecal][iptbin]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
      HoEUnsplitEle[iecal][iptbin]->Sumw2();
      if(iecal<2) {
        sprintf(histo,"sigmaIEtaIEtaUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
        sigmaIEtaIEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEBMin, sigmaIEtaIEtaEBMax);
        sprintf(histo,"sigmaIPhiIPhiUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
        sigmaIPhiIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEBMin, sigmaIPhiIPhiEBMax);
      } else {
        sprintf(histo,"sigmaIEtaIEtaUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
        sigmaIEtaIEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEEMin, sigmaIEtaIEtaEEMax);
        sprintf(histo,"sigmaIPhiIPhiUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
        sigmaIPhiIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEEMin, sigmaIPhiIPhiEEMax);
      }
      sigmaIEtaIEtaUnsplitEle[iecal][iptbin]->Sumw2();
      sigmaIPhiIPhiUnsplitEle[iecal][iptbin]->Sumw2();
      sprintf(histo,"fBremUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
      fBremUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, fBremMin, fBremMax);
      fBremUnsplitEle[iecal][iptbin]->Sumw2();
      sprintf(histo,"lhUnsplit_%s_subdet%d_ptbin%d",hypothesis,iecal,iptbin);
      lhUnsplitEle[iecal][iptbin]     = new TH1F(histo, histo, nbins, lhMin, lhMax);
      lhUnsplitEle[iecal][iptbin]->Sumw2();

      // iclass = 0: 0 - brem clusters
      // iclass = 1: >=1 - brem clusters
      for(int iclass=0; iclass<2; iclass++) {
	sprintf(histo,"dPhiClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
	dPhiClassEle[iecal][iptbin][iclass]     = new TH1F(histo, histo, nbins, dPhiMin, dPhiMax);
	dPhiClassEle[iecal][iptbin][iclass]->Sumw2();
	sprintf(histo,"dEtaClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
	dEtaClassEle[iecal][iptbin][iclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	dEtaClassEle[iecal][iptbin][iclass]->Sumw2();
	sprintf(histo,"EoPClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
	EoPClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, EoPMin, EoPMax);
	EoPClassEle[iecal][iptbin][iclass]->Sumw2();
	sprintf(histo,"OneOverEMinusOneOverPClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
	OneOverEMinusOneOverPClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, OneOverEMinusOneOverPMin, OneOverEMinusOneOverPMax);
	OneOverEMinusOneOverPClassEle[iecal][iptbin][iclass]->Sumw2();
	sprintf(histo,"HoEClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
	HoEClassEle[iecal][iptbin][iclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
	HoEClassEle[iecal][iptbin][iclass]->Sumw2();
        if(iecal<2) {
          sprintf(histo,"sigmaIEtaIEtaClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
          sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEBMin, sigmaIEtaIEtaEBMax);
          sprintf(histo,"sigmaIPhiIPhiClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
          sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEBMin, sigmaIPhiIPhiEBMax);
        } else {
          sprintf(histo,"sigmaIEtaIEtaClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
          sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaEEMin, sigmaIEtaIEtaEEMax);
          sprintf(histo,"sigmaIPhiIPhiClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
          sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiEEMin, sigmaIPhiIPhiEEMax);
        }
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass]->Sumw2();
	sigmaIPhiIPhiClassEle[iecal][iptbin][iclass]->Sumw2();
        sprintf(histo,"fBremClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
        fBremClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, fBremMin, fBremMax);
	fBremClassEle[iecal][iptbin][iclass]->Sumw2();
	sprintf(histo,"lhClass_%s_subdet%d_ptbin%d_class%d",hypothesis,iecal,iptbin,iclass);
	lhClassEle[iecal][iptbin][iclass]     = new TH1F(histo, histo, nbins, lhMin, lhMax);
	lhClassEle[iecal][iptbin][iclass]->Sumw2();
      }
    }
  }
}

float sPlotsPdfsComparison::likelihoodRatio(int isMc, ElectronLikelihood &lh) {
  LikelihoodMeasurements measurements;
  if(isMc) {
    measurements.pt = f_pt;
    if(fabs(f_eta)<1.0) measurements.subdet = 0;
    else if (fabs(f_eta)<1.479) measurements.subdet = 1;
    else measurements.subdet = 2;    
    measurements.deltaPhi = f_dphi;
    measurements.deltaEta = f_deta;
    measurements.eSeedClusterOverPout = -1;
    measurements.eSuperClusterOverP = f_eop;
    measurements.hadronicOverEm = f_hoe;
    measurements.sigmaIEtaIEta = f_see;
    measurements.sigmaIPhiIPhi = f_spp;
    measurements.fBrem = f_fbrem;
    measurements.nBremClusters = f_nbrem;
  } else {
    measurements.pt = pt;
    if(fabs(eta)<1.0) measurements.subdet = 0;
    else if (fabs(eta)<1.479) measurements.subdet = 1;
    else measurements.subdet = 2;    
    measurements.deltaPhi = dphi;
    measurements.deltaEta = deta;
    measurements.eSeedClusterOverPout = -1;
    measurements.eSuperClusterOverP = eop;
    measurements.hadronicOverEm = hoe;
    measurements.sigmaIEtaIEta = see;
    measurements.sigmaIPhiIPhi = spp;
    measurements.fBrem = fbrem;
    measurements.nBremClusters = nbrem;
  }
  return lh.result(measurements);
}
