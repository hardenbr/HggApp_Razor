#include <iostream>
#include <string> 
#include <math.h>

#include "TLorentzVector.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "EgammaAnalysisTools/include/HZZ4LElectronSelector.hh"
#include "EgammaAnalysisTools/include/eIDCiChzzSelector.hh"

#include "cajun/json/reader.h"
#include "cajun/json/elements.h"
#include <CommonTools/include/TriggerMask.hh>

using namespace bits;
using namespace std;

HZZ4LElectronSelector::HZZ4LElectronSelector(TTree *tree)
  : Egamma(tree) {

  // configuring the electron BDT
  fMVAHWW = new ElectronIDMVA();
  fMVAHWW->Initialize("BDTG method",
                      "elebdtweights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
                      "elebdtweights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml" ,                
                      ElectronIDMVA::kWithIPInfo);

  fMVAHWWWithIso = new ElectronIDMVA();
  fMVAHWWWithIso->Initialize("BDTG method",
			     "elebdtweights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml",
			     "elebdtweights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml" ,                
			     ElectronIDMVA::kIDIsoCombined);
  
  // configuring the electron BDT for H->ZZ
  fMVAHZZDanV0 = new ElectronIDMVAHZZ();
  fMVAHZZSiV0 = new ElectronIDMVAHZZ();
  fMVAHZZSiV1 = new ElectronIDMVAHZZ();
  fMVAHZZSiDanV2 = new ElectronIDMVAHZZ();

  // New H->ZZ unbiased DATA training, Daniele's variables
  fMVAHZZDanV0->Initialize("BDTSimpleCat",
                           "elebdtweights/DanieleMVA_BDTCat_BDTG_DanV0.weights.xml",
                           ElectronIDMVAHZZ::kBDTDanV0);

  // New H->ZZ unbiased DATA training, Si's HWW 2011 variables
  fMVAHZZSiV0->Initialize("BDTSimpleCat",
                          "elebdtweights/DanieleMVA_BDTCat_BDTG_SiV0.weights.xml",
                          ElectronIDMVAHZZ::kBDTSiV0);

  // New H->ZZ unbiased DATA training, Si's HWW 2012 variables
  fMVAHZZSiV1->Initialize("BDTSimpleCat",
                          "elebdtweights/DanieleMVA_BDTCat_BDTG_SiV1.weights.xml",
                          ElectronIDMVAHZZ::kBDTSiV1);

  // New H->ZZ unbiased DATA training, Daniele's + Si's variables 
  fMVAHZZSiDanV2->Initialize("BDTSimpleCat",
                             "elebdtweights/DanieleMVA_BDTCat_BDTG_SiDanV2.weights.xml",
                             ElectronIDMVAHZZ::kBDTSiDanV2);

  // configuring the electron BDT for H->WW
  fMVAHWWDanV0 = new ElectronIDMVAHZZ();
  fMVAHWWSiV0 = new ElectronIDMVAHZZ();
  fMVAHWWSiV1 = new ElectronIDMVAHZZ();
  fMVAHWWSiDanV2 = new ElectronIDMVAHZZ();

  // New H->WW unbiased DATA training, Daniele's variables
  fMVAHWWDanV0->Initialize("BDTSimpleCat",
                           "elebdtweights/DanieleMVA_DenomHWW_BDTCat_BDTG_DanV0.weights.xml",
                           ElectronIDMVAHZZ::kBDTHWWDanV0);

  // New H->WW unbiased DATA training, Si's HWW 2011 variables
  fMVAHWWSiV0->Initialize("BDTSimpleCat",
                          "elebdtweights/DanieleMVA_DenomHWW_BDTCat_BDTG_SiV0.weights.xml",
                          ElectronIDMVAHZZ::kBDTHWWSiV0);

  // New H->WW unbiased DATA training, Si's HWW 2012 variables
  fMVAHWWSiV1->Initialize("BDTSimpleCat",
                          "elebdtweights/DanieleMVA_DenomHWW_BDTCat_BDTG_SiV1.weights.xml",
                          ElectronIDMVAHZZ::kBDTHWWSiV1);

  // New H->WW unbiased DATA training, Daniele's + Si's variables 
  fMVAHWWSiDanV2->Initialize("BDTSimpleCat",
                             "elebdtweights/DanieleMVA_DenomHWW_BDTCat_BDTG_SiDanV2.weights.xml",
                             ElectronIDMVAHZZ::kBDTHWWSiDanV2);

}

HZZ4LElectronSelector::~HZZ4LElectronSelector() { 
  
}


// PDF for probe electrons within the acceptance and loose isolated in the tracker
// the tag is the one with the best match to the Z
// the tag must be within the acceptance, tracker isolated and loose identified
void HZZ4LElectronSelector::Loop(const char *treefilesuffix) {
  
  if(fChain == 0) return;

  Utils anaUtils;
  
  char treename[200];
  sprintf(treename,"%s.root",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addElectronIdBits();
  reducedTree.addIsolations();
  reducedTree.addRunInfos();
  reducedTree.addMore();
  reducedTree.addTrackMomenta();

  // counters
  int allevents   = 0;
  int etaprobes   = 0;
  int ptprobes    = 0;
  int matchprobes = 0;
  int goodevents  = 0;

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();

  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    allevents++;
    
    electrons.clear();
    for(int iele=0; iele<nEle; iele++) {
      TLorentzVector electron(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);

      if( fabs(etaEle[iele])>2.5 ) continue;
      etaprobes++;
      
      if( electron.Pt()<5 ) continue;
      ptprobes++;

      if( !mcMatches(iele) ) continue;
      matchprobes++;

      electrons.push_back(iele);
    }

    if(electrons.size()==0) continue;
    goodevents++;

    for(int i=0; i<(int)electrons.size(); ++i) {
      
      int probe = electrons[i];
      TLorentzVector probeP4(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);

      // various about probe
      int charge = chargeEle[probe];
      float pt   = probeP4.Pt();
      float eta  = etaEle[probe];
      float phi  = phiEle[probe];
      
      // some eleID variables
      float HoE, s1s9, s9s25, phiwidth, etawidth, deta, dphi, fbrem, see, spp, eleopout, eopout, eop, eseedopin, nbrems, recoFlag, EleSCEta, EleSCPhi;
      float EleSCEt,EtaSeed,PhiSeed,ESeed,IEtaSeed,IPhiSeed,EtaCrySeed,PhiCrySeed,IEtaCrySeed,IPhiCrySeed;
      float oneoveremoneoverp, eledeta, d0, ip3d, ip3ds, kfhits, kflayers, kfchi2, e1x5e5x5, dcot, dist;
      float detacalo, dphicalo, sep, dz, gsfchi2, emax, etop, ebottom, eleft, eright,
	e2nd, e2x5right, e2x5left, e2x5top, e2x5bottom, 
	e2x5max, e1x5, e2x2, e3x3, e5x5, r9, nclu,
	scenergy, scrawenergy, scesenergy;
      
      int gsfTrack = gsfTrackIndexEle[probe];   
      int kfTrack = trackIndexEle[probe];

      // different p estimations
      float pcomb=probeP4.Vect().Mag();
      TVector3 p3ModeGsf(pxModeGsfTrack[gsfTrack],pyModeGsfTrack[gsfTrack],pzModeGsfTrack[gsfTrack]);
      float pmodegsf=p3ModeGsf.Mag();
      TVector3 p3MeanGsf(pxGsfTrack[gsfTrack],pyGsfTrack[gsfTrack],pzGsfTrack[gsfTrack]);
      float pmeangsf=p3MeanGsf.Mag();
      TVector3 p3MeanKf(pxTrack[kfTrack],pyTrack[kfTrack],pzTrack[kfTrack]);
      float pmeankf=p3MeanKf.Mag();
      float pterrorgsf = ptErrorGsfTrack[gsfTrack];
      float pterrorkf = (kfTrack>-1) ? ptErrorTrack[kfTrack] : -1.0;
      float perrorele = trackMomentumErrorEle[probe];

      double gsfsign   = (-eleDxyPV(probe,0) >=0 ) ? 1. : -1.;
      bool matchConv = hasMatchedConversionEle[probe];
      
      d0 = gsfsign * transvImpactParGsfTrack[gsfTrack];
      dz = eleDzPV(probe,0);
      ip3d = gsfsign * impactPar3DGsfTrack[gsfTrack];
      ip3ds = ip3d/impactPar3DErrorGsfTrack[gsfTrack];
      kfchi2 = (kfTrack>-1) ? trackNormalizedChi2Track[kfTrack] : 0.0;
      kflayers = (kfTrack>-1) ? trackerLayersWithMeasurementTrack[kfTrack] : -1.0;
      kfhits = (kfTrack>-1) ? trackValidHitsTrack[kfTrack] : -1.0;
      gsfchi2 = trackNormalizedChi2GsfTrack[gsfTrack];
      int misshits = expInnerLayersGsfTrack[gsfTrack];
      dcot = convDistEle[probe];
      dist = convDcotEle[probe];      
      bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[probe], isEcalDriven);
      int ecalseed;
      HoE = hOverEEle[probe];
      eledeta = deltaEtaEleClusterTrackAtCaloEle[probe];
      deta = deltaEtaAtVtxEle[probe];
      dphi = deltaPhiAtVtxEle[probe];
      detacalo = deltaEtaAtCaloEle[probe];
      dphicalo = deltaPhiAtCaloEle[probe];
      fbrem = fbremEle[probe];
      nbrems = nbremsEle[probe];
      eleopout = eEleClusterOverPoutEle[probe];
      eopout = eSeedOverPoutEle[probe];
      eop = eSuperClusterOverPEle[probe];
      if(ecaldriven) {
	ecalseed = 1;
	int sc = superClusterIndexEle[probe];
	float seedEnergy = seedClusterEnergySC[sc];
        nclu = float(nBCSC[sc]);
	s1s9 = eMaxSC[sc]/eMaxSC[sc];
	s9s25 = e3x3SC[sc]/e5x5SC[sc];
	e1x5e5x5 = (e5x5SC[sc] - e1x5SC[sc])/e5x5SC[sc];
	phiwidth = phiWidthSC[sc];
	etawidth = etaWidthSC[sc];
	see = sqrt(covIEtaIEtaSC[sc]);
	sep = covIEtaIPhiSC[sc]/(sqrt(covIEtaIEtaSC[sc])*sqrt(covIPhiIPhiSC[sc]));
	spp = sqrt(covIPhiIPhiSC[sc]);
	oneoveremoneoverp = 1./energySC[sc]  - 1./probeP4.Vect().Mag();
	eseedopin = seedEnergy/p3ModeGsf.Mag();
	emax = eMaxSC[sc];
	etop = eTopSC[sc];
	ebottom = eBottomSC[sc];
	eleft = eLeftSC[sc];
	eright = eRightSC[sc];
	e2nd = e2ndSC[sc];
	e2x5right = e2x5RightSC[sc];
	e2x5left = e2x5LeftSC[sc];
	e2x5top = e2x5TopSC[sc];
	e2x5bottom = e2x5BottomSC[sc];
	e2x5max = e2x5MaxSC[sc];
	e1x5 = e1x5SC[sc];
	e2x2 = e2x2SC[sc];
	e3x3 = e3x3SC[sc];
	e5x5 = e5x5SC[sc];
	r9 = e3x3SC[sc]/rawEnergySC[sc];
	recoFlag = recoFlagSC[sc];
	EleSCEta = etaSC[sc];
	EleSCPhi = phiSC[sc];            
        EleSCEt = energySC[sc]*fabs(sin(thetaSC[sc]));
	scenergy = energySC[sc];
	scrawenergy = rawEnergySC[sc];
	scesenergy = esEnergySC[sc];
        int seedclu = indexSeedBC(sc,ecaldriven);
        EtaSeed=etaBC[seedclu];
        PhiSeed=phiBC[seedclu];
        ESeed=energyBC[seedclu];
        IEtaSeed=iEtaBC[seedclu];
        IPhiSeed=iPhiBC[seedclu];
        EtaCrySeed=etaCrystalBC[seedclu];
        PhiCrySeed=phiCrystalBC[seedclu];
        IEtaCrySeed=seedXSC[sc];
        IPhiCrySeed=seedYSC[sc];
      } else {
	ecalseed = 0;
	int sc = PFsuperClusterIndexEle[probe];
	if(sc>-1) {
	  float seedEnergy = seedClusterEnergyPFSC[sc];
          nclu = float(nBCPFSC[sc]);
	  s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
	  s1s9 = eMaxPFSC[sc]/eMaxPFSC[sc];
	  e1x5e5x5 = (e5x5PFSC[sc] - e1x5PFSC[sc])/e5x5PFSC[sc];
	  phiwidth = phiWidthPFSC[sc];
	  etawidth = etaWidthPFSC[sc];
	  see = sqrt(covIEtaIEtaPFSC[sc]);
	  sep = covIEtaIPhiPFSC[sc]/(sqrt(covIEtaIEtaPFSC[sc])*sqrt(covIPhiIPhiPFSC[sc]));
	  spp = sqrt(covIPhiIPhiPFSC[sc]);
	  oneoveremoneoverp = 1./energyPFSC[sc]  - 1./probeP4.Vect().Mag();
	  eseedopin = seedEnergy/p3ModeGsf.Mag();
	  emax = eMaxPFSC[sc];
	  etop = eTopPFSC[sc];
	  ebottom = eBottomPFSC[sc];
	  eleft = eLeftPFSC[sc];
	  eright = eRightPFSC[sc];
	  e2nd = e2ndPFSC[sc];
	  e2x5right = e2x5RightPFSC[sc];
	  e2x5left = e2x5LeftPFSC[sc];
	  e2x5top = e2x5TopPFSC[sc];
	  e2x5bottom = e2x5BottomPFSC[sc];
	  e2x5max = e2x5MaxPFSC[sc];
	  e1x5 = e1x5PFSC[sc];
	  e2x2 = e2x2PFSC[sc];
	  e3x3 = e3x3PFSC[sc];
	  e5x5 = e5x5PFSC[sc];
	  r9 = e3x3PFSC[sc]/rawEnergyPFSC[sc];
	  recoFlag = recoFlagPFSC[sc];
	  EleSCEta = etaPFSC[sc];
	  EleSCPhi = phiPFSC[sc];     
          EleSCEt = energyPFSC[sc]*fabs(sin(thetaPFSC[sc]));
	  scenergy = energyPFSC[sc];
	  scrawenergy = rawEnergyPFSC[sc];
	  scesenergy = esEnergyPFSC[sc];
          int seedclu = indexSeedBC(sc,ecaldriven);
          EtaSeed=etaPFBC[seedclu];
          PhiSeed=phiPFBC[seedclu];
          ESeed=energyPFBC[seedclu];
          IEtaSeed=iEtaPFBC[seedclu];
          IPhiSeed=iPhiPFBC[seedclu];
          EtaCrySeed=etaCrystalPFBC[seedclu];
          PhiCrySeed=phiCrystalPFBC[seedclu];
          IEtaCrySeed=seedXPFSC[sc];
          IPhiCrySeed=seedYPFSC[sc];
	} else {
	  s9s25 = 999.;
	  see = 999.;
	  spp = 999.;
	}
      }

      // this only happens in 42X EE because of a bug in the EcalClusterTools
      if(isnan(see) || isnan(spp) || isnan(sep)) continue;

      // CiC...
      eIDCiChzzSelector cicsel;
      int cic[5];
      int iecal = (fabs(eta)<1.479) ? 0 : 1;
      for(int i=0; i<5; i++) cic[i] = cicsel.ElectronId_V03(pt,EleSCEta,see,eop,eseedopin,fbrem,
							    dr03TkSumPtEle[probe],
							    dr03EcalRecHitSumEtEle[probe] - rhoFastjet * Aeff_ecal_dr03[iecal],
							    dr03HcalTowerSumEtEle[probe] - rhoFastjet * Aeff_hcal_dr03[iecal],
							    d0,misshits,deta,dphi,HoE,dcot,
							    dist,!ecalseed,i,false);

      // some MVAs...
      float pfmva = pflowMVAEle[probe];
      float lh=eleIdLikelihoodEle[probe];
      float hwwbdts[2];
      hwwbdts[0] = eleBDT(fMVAHWW,probe);
      hwwbdts[1] = eleBDTWithIso(fMVAHWWWithIso,probe);
      float hzzbdts[4];
      hzzbdts[0] = eleBDT(fMVAHZZDanV0,probe);
      hzzbdts[1] = eleBDT(fMVAHZZSiV0,probe);
      hzzbdts[2] = eleBDT(fMVAHZZSiV1,probe);
      // hzzbdts[3] = eleBDT(fMVAHZZSiDanV2,probe);
      hzzbdts[3] = mvaidnontrigEle[probe];
      float newhwwbdts[4];
      newhwwbdts[0] = eleBDT(fMVAHWWDanV0,probe);
      newhwwbdts[1] = eleBDT(fMVAHWWSiV0,probe);
      newhwwbdts[2] = eleBDT(fMVAHWWSiV1,probe);
      // newhwwbdts[3] = eleBDT(fMVAHWWSiDanV2,probe);
      newhwwbdts[3] = mvaidtrigEle[probe];

      // isolations
      float chaPfIso[8], phoPfIso[8], neuPfIso[8];
      chaPfIso[0]=pfCandChargedIso01Ele[probe];
      phoPfIso[0]=pfCandPhotonIso01Ele[probe];
      neuPfIso[0]=pfCandNeutralIso01Ele[probe];
      
      chaPfIso[1]=pfCandChargedIso02Ele[probe];
      phoPfIso[1]=pfCandPhotonIso02Ele[probe];
      neuPfIso[1]=pfCandNeutralIso02Ele[probe];
      
      chaPfIso[2]=pfCandChargedIso03Ele[probe];
      phoPfIso[2]=pfCandPhotonIso03Ele[probe];
      neuPfIso[2]=pfCandNeutralIso03Ele[probe];

      chaPfIso[3]=pfCandChargedIso04Ele[probe];
      phoPfIso[3]=pfCandPhotonIso04Ele[probe];
      neuPfIso[3]=pfCandNeutralIso04Ele[probe];

      chaPfIso[4]=pfCandChargedIso05Ele[probe];
      phoPfIso[4]=pfCandPhotonIso05Ele[probe];
      neuPfIso[4]=pfCandNeutralIso05Ele[probe];
      
      chaPfIso[5]=pfCandChargedIso06Ele[probe];
      phoPfIso[5]=pfCandPhotonIso06Ele[probe];
      neuPfIso[5]=pfCandNeutralIso06Ele[probe];
      
      chaPfIso[6]=pfCandChargedIso07Ele[probe];
      phoPfIso[6]=pfCandPhotonIso07Ele[probe];
      neuPfIso[6]=pfCandNeutralIso07Ele[probe];
      
      chaPfIso[7]=pfCandChargedDirIso04Ele[probe];
      phoPfIso[7]=pfCandPhotonDirIso04Ele[probe];
      neuPfIso[7]=pfCandNeutralDirIso04Ele[probe];
      
      float trkIso[2], ecalIso[2], hcalIso[2];
      trkIso[0]=dr03TkSumPtEle[probe];
      trkIso[1]=dr04TkSumPtEle[probe];
      ecalIso[0]=dr03EcalRecHitSumEtEle[probe];
      ecalIso[1]=dr04EcalRecHitSumEtEle[probe];
      hcalIso[0]=dr03HcalTowerSumEtFullConeEle[probe];
      hcalIso[1]=dr04HcalTowerSumEtFullConeEle[probe];
      
      bool isEleEB= anaUtils.fiducialFlagECAL(fiducialFlagsEle[probe], isEB);
      bool isEleEE= anaUtils.fiducialFlagECAL(fiducialFlagsEle[probe], isEE);

      // fill the reduced tree
      reducedTree.fillVariables(eleopout,eopout,eop,HoE,deta,dphi,s9s25,s1s9,see,spp,fbrem,
				nbrems,misshits,dcot,dist,pt,eta,charge,phiwidth,etawidth,
				oneoveremoneoverp,eledeta,d0,ip3d,ip3ds,kfhits,kflayers,kfchi2,e1x5e5x5,ecalseed,matchConv,
                                isEleEB,isEleEE);
      reducedTree.fillVariables2(detacalo, dphicalo, sep, dz, gsfchi2, emax, etop, ebottom, eleft, eright,
				 e2nd, e2x5right, e2x5left, e2x5top, e2x5bottom, 
				 e2x5max, e1x5, e2x2, e3x3, e5x5, r9, nclu,
				 phi, scenergy, scrawenergy, scesenergy, eseedopin);
      reducedTree.fillIsolations(trkIso,ecalIso,hcalIso, 
				 pfCombinedIsoEle[probe],
				 chaPfIso, neuPfIso, phoPfIso);
      reducedTree.fillCluterInfos(EleSCEt,EleSCEta,EleSCPhi,EtaSeed,PhiSeed,ESeed,IEtaSeed,IPhiSeed,EtaCrySeed,PhiCrySeed,IEtaCrySeed,IPhiCrySeed);
      reducedTree.fillMore(nPV,rhoFastjet,hwwbdts,newhwwbdts,hzzbdts,pfmva,lh);
      reducedTree.fillTrackMomenta(pcomb,pmodegsf,pmeangsf,pmeankf,pterrorgsf,pterrorkf,perrorele);
      reducedTree.fillFakeRateDenomBits(-1.,isDenomFake(probe),isDenomFake_smurfs(probe));
      reducedTree.fillBDTBasedIDBits(passEleBDT(pt,EleSCEta,hwwbdts[0]));
      reducedTree.fillCiCBasedIDBits(cic);
      reducedTree.fillRunInfos(runNumber, lumiBlock, eventNumber, nPU, 1);
      reducedTree.store();
    } // possible probes
    
  } // loop over events
  
  cout << "statistics from Tag and Probe: " << endl;
  cout << "allevents   = " << allevents << endl;
  cout << "etaprobes   = " << etaprobes << endl;
  cout << "ptprobes    = " << ptprobes << endl;
  cout << "matchprobes = " << matchprobes << endl;
  cout << "goodevents  = " << goodevents << endl;

  reducedTree.save();
}


// denominator for fake rate: for HtoWW, egamma triggers
int HZZ4LElectronSelector::isDenomFake(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);

  // match with the HLT firing candidates
  //  bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
  //  if (!HLTmatch) isGoodDenom = false;
  
  // acceptance for the fake electron
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )        isGoodDenom = false;

  // barrel or endcap
  bool isEleEB = anaUtils.fiducialFlagECAL(fiducialFlagsEle[theEle], isEB);    

  // taking shower shape                                                                                                             
  int sc;
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  float thisSigmaIeIe = -1.;
  if (ecalDriven) {
    sc = superClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaSC[sc]);
  }
  if (!ecalDriven) {
    sc = PFsuperClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaPFSC[sc]);
  }
  if ( sc<0 ) { isGoodDenom = false; }

  // sigmaIetaIeta                                                                                                                   
  if ( isEleEB && thisSigmaIeIe>0.01) { isGoodDenom = false; }
  if (!isEleEB && thisSigmaIeIe>0.03) { isGoodDenom = false; }

  // H/E
  if ( isEleEB && hOverEEle[theEle]>0.12) isGoodDenom = false;   
  if (!isEleEB && hOverEEle[theEle]>0.10) isGoodDenom = false;
  
  // deltaEta
  if ( isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.007) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.009) ) isGoodDenom = false;
  
  // deltaPhi
  if ( isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.15) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.10) ) isGoodDenom = false;

  // isolation 
  float ecalIsol    = (dr03EcalRecHitSumEtEle[theEle])/p3Ele.Pt();
  float hcalIsol    = (dr03HcalTowerSumEtEle[theEle])/p3Ele.Pt();
  float trackerIsol = (dr03TkSumPtEle[theEle])/p3Ele.Pt();                
  if(ecalIsol>0.2)    isGoodDenom = false;
  if(hcalIsol>0.2)    isGoodDenom = false;
  if(trackerIsol>0.2) isGoodDenom = false;                                
   
  if(isGoodDenom) return 1;
  return 0;
}

// denominator for fake rate: for HtoWW, egamma triggers, same as smurfs
int HZZ4LElectronSelector::isDenomFake_smurfs(int theEle) {
  
  Utils anaUtils;
  bool isGoodDenom = true;
  
  TVector3 p3Ele(pxEle[theEle], pyEle[theEle], pzEle[theEle]);
  
  // match with the HLT firing candidates
  // bool HLTmatch = triggerMatch(p3Ele.Eta(),p3Ele.Phi(),0.2);
  // if (!HLTmatch) isGoodDenom = false;
  
  // acceptance for the fake electron
  if( fabs(p3Ele.Eta()) > 2.5 ) isGoodDenom = false;
  if( p3Ele.Pt() < 10. )        isGoodDenom = false;

  // taking shower shape                                                                                                             
  int sc;
  bool ecalDriven = anaUtils.electronRecoType(recoFlagsEle[theEle], bits::isEcalDriven);
  float thisSigmaIeIe = -1.;
  float scEta = -1.;                                                               
  if (ecalDriven) {
    sc = superClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaSC[sc]);
    scEta = etaSC[sc];
  }
  if (!ecalDriven) {
    sc = PFsuperClusterIndexEle[theEle];
    thisSigmaIeIe = sqrt(covIEtaIEtaPFSC[sc]);
    scEta = etaPFSC[sc];
  }
  if ( sc<0 ) { isGoodDenom = false; }

  // barrel or endcap
  bool isEleEB = false;
  if (fabs(scEta)<1.479) isEleEB = true;   
  
  // sigmaIetaIeta                                                                                                                   
  if ( isEleEB && thisSigmaIeIe>0.01) { isGoodDenom = false; }
  if (!isEleEB && thisSigmaIeIe>0.03) { isGoodDenom = false; }

  // isolation
  float ecalIsolAbs = 0.0;
  if ( isEleEB ) ecalIsolAbs = max(0.0,dr03EcalRecHitSumEtEle[theEle]-1.0);
  else ecalIsolAbs = dr03EcalRecHitSumEtEle[theEle];
  float ecalIsol = ecalIsolAbs/p3Ele.Pt(); 
  float hcalIsol    = (dr03HcalTowerSumEtEle[theEle])/p3Ele.Pt();
  float trackerIsol = (dr03TkSumPtEle[theEle])/p3Ele.Pt();                
  if(ecalIsol>0.2)    isGoodDenom = false;
  if(hcalIsol>0.2)    isGoodDenom = false;
  if(trackerIsol>0.2) isGoodDenom = false;                                

  // H/E
  if ( isEleEB && hOverEEle[theEle]>0.12) isGoodDenom = false;   
  if (!isEleEB && hOverEEle[theEle]>0.10) isGoodDenom = false;
  
  // deltaEta
  if ( isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.007) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaEtaAtVtxEle[theEle])>0.009) ) isGoodDenom = false;
  
  // deltaPhi
  if ( isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.15) ) isGoodDenom = false;
  if (!isEleEB && (fabs(deltaPhiAtVtxEle[theEle])>0.10) ) isGoodDenom = false;
  
  // full conversion rejection 
  int gsf = gsfTrackIndexEle[theEle];
  int missHits = expInnerLayersGsfTrack[gsf];
  bool matchConv = hasMatchedConversionEle[theEle];
  if (missHits>0 || matchConv) isGoodDenom = false;

  // impact parameter cuts 
  float dxyEle = transvImpactParGsfTrack[gsf];
  float dzEle  = PVzPV[0] - trackVzGsfTrack[gsf];
  if (fabs(dxyEle)>0.02) isGoodDenom = false;
  if (fabs(dzEle)>0.10)  isGoodDenom = false;

  if(isGoodDenom) return 1;
  return 0;
}


bool HZZ4LElectronSelector::mcMatches(int probe) {
  
  bool probematch=true;
  bool tagmatch=true;

  std::vector<int> ep,em;
  
  for(int imc=0; imc<20; ++imc) {
    if(idMc[imc]==11 && abs(idMc[mothMc[imc]])>=23 && abs(idMc[mothMc[imc]])<=24 && abs(idMc[mothMc[mothMc[imc]]])==25 && em.size()<2) em.push_back(imc);
    else if(idMc[imc]==-11 && abs(idMc[mothMc[imc]])>=23 && abs(idMc[mothMc[imc]])<=24 && abs(idMc[mothMc[mothMc[imc]]])==25 && ep.size()<2) ep.push_back(imc);
  }
  
  TVector3 probeP(pxEle[probe],pyEle[probe],pzEle[probe]);
  bool matches=false;
  if(chargeEle[probe]>0) {
    for(int imc=0;imc<(int)ep.size();++imc) {
      int probemc=ep[imc];
      TVector3 probemcP;
      probemcP.SetMagThetaPhi(pMc[probemc],thetaMc[probemc],phiMc[probemc]);
      if(probemcP.DeltaR(probeP)<0.1) matches=true;
    }
  } else {
    for(int imc=0;imc<(int)em.size();++imc) {
      int probemc=em[imc];
      TVector3 probemcP;
      probemcP.SetMagThetaPhi(pMc[probemc],thetaMc[probemc],phiMc[probemc]);
      if(probemcP.DeltaR(probeP)<0.1) matches=true;
    }
  }
  return matches;
}
