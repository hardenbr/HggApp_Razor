#include <HggPhotonID.hh>
#include "ReadConfig.hh"

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
using namespace std;
using namespace TMVA;

#define debugPhotonID 0

HggPhotonID::HggPhotonID():
  doECALIso(true)
{
}

void HggPhotonID::Init(){

  ReadConfig cfg;
  if(cfg.read(configFile)!=0){
    cout << "ERROR: Could not read PhotonID Config!";
    valid = false;
    return;
  }

  /*  weightFile_IdEB_2011 = cfg.getParameter("weightFile_IdEB_2011");
  weightFile_IdEE_2011 = cfg.getParameter("weightFile_IdEE_2011");
  weightFile_IdEB_2012 = cfg.getParameter("weightFile_IdEB_2012");
  weightFile_IdEE_2012 = cfg.getParameter("weightFile_IdEE_2012");
  */ //RAZOR UNNCESSARY
  methodName_Id  = cfg.getParameter("methodName_Id");

  version = cfg.getParameter("PhotonIDVersion");

  if(cfg.getParameter("doEcalIso").compare("no")==0) doECALIso=false;
  else doECALIso=true;

  cout << "doECALIso: " << doECALIso <<endl;

  photonMVA_EB_2011 = new TMVA::Reader( "!Color:!Silent" );;
  photonMVA_EE_2011 = new TMVA::Reader( "!Color:!Silent" );;
  photonMVA_EB_2012 = new TMVA::Reader( "!Color:!Silent" );;
  photonMVA_EE_2012 = new TMVA::Reader( "!Color:!Silent" );;

  InputDists["pfPhotonIsooet"] = new TH1F("pfPhotonIsooet","",200,0,100);
  InputDists["pfChargedIsooet"] = new TH1F("pfChargedIsooet","",200,0,100);
  InputDists["pfChargedIsoWorstoet"] = new TH1F("pfChargedIsoWorstoet","",200,0,100);
  InputDists["isosumoet"] = new TH1F("isosumoet","",200,0,100);
  InputDists["isosumoetbad"] = new TH1F("isosumoetbad","",200,0,100);
  InputDists["isosumoetPF"] = new TH1F("isosumoetPF","",200,0,100);
  InputDists["isosumoetbadPF"] = new TH1F("isosumoetbadPF","",200,0,100);

  //  this->setupTMVA();//not necessary for RAZOR
}

void HggPhotonID::setVertices(int nPV, float* xPV, float *yPV, float *zPV){
  vertices.clear();
  for(int i=0;i<nPV;i++){
    vertices.push_back(TVector3(xPV[i],yPV[i],zPV[i]));
  }
}

void HggPhotonID::fillVariables(VecbosPho* pho, int nVertex, float rhoFastJet,int selVtxIndex){
  selVtxPos = vertices.at(selVtxIndex);

  eT = pho->p4FromVtx(selVtxPos,pho->finalEnergy,false).Et();
  if(debugPhotonID) cout << "eT: " << eT << endl;
  hoe = pho->HoverE;
  sigietaieta=pho->SC.sigmaIEtaIEta;
  r9 = pho->SC.r9;
  ecalisodr03 = pho->dr03EcalRecHitSumEtCone;
  ecalisodr04 = pho->dr04EcalRecHitSumEtCone;
  hcalisodr03 = pho->dr03HcalTowerSumEtCone;
  hcalisodr04 = pho->dr04HcalTowerSumEtCone;
  nVertexf=nVertex;
  etasc = pho->SC.eta;
  scetawidth = pho->SC.etaWidth;
  scphiwidth = pho->SC.phiWidth;
                                                                                                                             
  pfPhotonIso03 = pho->dr03PhotonPFIso;
  pfPhotonIso03oet = pfPhotonIso03*50./eT;
  pfPhotonIso04 = pho->dr04PhotonPFIso;
  pfPhotonIso04oet = pfPhotonIso04*50./eT;
  sigietaiphi = pho->SC.sigmaIEtaIPhi;
  s4Ratio = pho->SC.e2x2/pho->SC.e5x5;
  rho = rhoFastJet;                                                                                                                 
  sigRR = pho->SC.esEffSigRR;                

  //We need some computation for these:

  pfChargedIsoGood03 = pho->dr03ChargedHadronPFIso[selVtxIndex];
  float maxIso=0;
  int worstVtx=0;
  for(int i=0; i<pho->nPV;i++) 
    if(pho->dr04ChargedHadronPFIso[i] > maxIso) 
      { 
	maxIso = pho->dr04ChargedHadronPFIso[i];
	worstVtx=i;
      }
  pfChargedIsoBad04  = maxIso;
  eTBad = pho->p4FromVtx(vertices.at(worstVtx),pho->finalEnergy,false).Et();

  pfChargedIsoGood03oet = pfChargedIsoGood03*50./eT;
  pfChargedIsoBad04oet  = pfChargedIsoBad04*50./eT;  
  
  trkisooet = pho->dr03TrackIso[selVtxIndex]*50./eT;
  isosum = (pho->dr03TrackIso[selVtxIndex]
	    + pho->dr03EcalRecHitSumEtCone 
	    + pho->dr04HcalTowerSumEtCone + isoSumConst - rho*rhoFac);
  isosumoet = isosum*50./eT;
  isosumbad = (*std::max_element(pho->dr04TrackIso,pho->dr04TrackIso+pho->nPV)
	       + pho->dr03EcalRecHitSumEtCone 
	       + pho->dr04HcalTowerSumEtCone + isoSumConst - rho*rhoFac);
  isosumoetbad = isosumbad*50./eT; 
  isosumPF = (pfChargedIsoGood03 + pfPhotonIso03 + isoSumConstPF - rho*rhoFac); 
  isosumbadPF = (pfChargedIsoBad04 + pfPhotonIso04 + isoSumConstPF - rho*rhoFacBad); 
  isosumoetPF = isosumPF*50./eT; 
  isosumoetbadPF = isosumbadPF*50./eTBad;

  InputDists["pfPhotonIsooet"]->Fill(pfPhotonIso03oet);
  InputDists["pfChargedIsooet"]->Fill(pfChargedIsoGood03oet);
  InputDists["pfChargedIsoWorstoet"]->Fill(pfChargedIsoBad04oet);
  InputDists["isosumoet"]->Fill(isosumoet);
  InputDists["isosumoetbad"]->Fill(isosumoetbad);
  InputDists["isosumoetPF"]->Fill(isosumoetPF);
  InputDists["isosumoetbadPF"]->Fill(isosumoetbadPF);

}

void HggPhotonID::setupTMVA(){

  //
  //   2011 Photon ID MVA
  //

  //setup the ID MVAs:
  ///EB
  photonMVA_EB_2011->AddVariable("HoE",&hoe);         
  photonMVA_EB_2011->AddVariable("covIEtaIEta",&sigietaieta);         
  photonMVA_EB_2011->AddVariable("tIso1abs",&isosumoet);              
  photonMVA_EB_2011->AddVariable("tIso3abs",&trkisooet);              
  photonMVA_EB_2011->AddVariable("tIso2abs",&isosumoetbad);           
  photonMVA_EB_2011->AddVariable("R9",&r9);           
  photonMVA_EB_2011->AddVariable("absIsoEcal",&ecalisodr03);          
  photonMVA_EB_2011->AddVariable("absIsoHcal",&hcalisodr04);          
  photonMVA_EB_2011->AddVariable("NVertexes",&nVertexf);              
  photonMVA_EB_2011->AddVariable("ScEta",&etasc);     
  photonMVA_EB_2011->AddVariable("EtaWidth",&scetawidth);             
  photonMVA_EB_2011->AddVariable("PhiWidth",&scphiwidth);             
       
  ///EE
  photonMVA_EE_2011->AddVariable("HoE",&hoe);         
  photonMVA_EE_2011->AddVariable("covIEtaIEta",&sigietaieta);         
  photonMVA_EE_2011->AddVariable("tIso1abs",&isosumoet);              
  photonMVA_EE_2011->AddVariable("tIso3abs",&trkisooet);              
  photonMVA_EE_2011->AddVariable("tIso2abs",&isosumoetbad);           
  photonMVA_EE_2011->AddVariable("R9",&r9);           
  photonMVA_EE_2011->AddVariable("absIsoEcal",&ecalisodr03);          
  photonMVA_EE_2011->AddVariable("absIsoHcal",&hcalisodr04);          
  photonMVA_EE_2011->AddVariable("NVertexes",&nVertexf);              
  photonMVA_EE_2011->AddVariable("ScEta",&etasc);     
  photonMVA_EE_2011->AddVariable("EtaWidth",&scetawidth);             
  photonMVA_EE_2011->AddVariable("PhiWidth",&scphiwidth);             
  
  //book MVAs:
  /*  photonMVA_EB_2011->BookMVA( methodName_Id, weightFile_IdEB_2011);
      photonMVA_EE_2011->BookMVA( methodName_Id, weightFile_IdEE_2011);*/ //RAZOR

  //
  //   2012 Photon ID MVA
  //

  photonMVA_EB_2012->AddVariable("myphoton_pfchargedisogood03",   &pfChargedIsoGood03oet );
  photonMVA_EB_2012->AddVariable("myphoton_pfchargedisobad03",   &pfChargedIsoBad04oet );
  photonMVA_EB_2012->AddVariable("myphoton_pfphotoniso03",   &pfPhotonIso03oet );

  photonMVA_EB_2012->AddVariable("myphoton_sieie",   &sigietaieta );
  photonMVA_EB_2012->AddVariable("myphoton_sieip",   &sigietaiphi );
  photonMVA_EB_2012->AddVariable("myphoton_etawidth",   &scetawidth );
  photonMVA_EB_2012->AddVariable("myphoton_phiwidth",   &scphiwidth );
  photonMVA_EB_2012->AddVariable("myphoton_r9",   &r9 );
  
  photonMVA_EB_2012->AddVariable("myphoton_s4ratio",   &s4Ratio );
  
  photonMVA_EB_2012->AddVariable("myphoton_SCeta",   &etasc );
  photonMVA_EB_2012->AddVariable("event_rho",   &rho );
  
  photonMVA_EE_2012->AddVariable("myphoton_pfchargedisogood03",   &pfChargedIsoGood03oet );
  photonMVA_EE_2012->AddVariable("myphoton_pfchargedisobad03",   &pfChargedIsoBad04oet );
  photonMVA_EE_2012->AddVariable("myphoton_pfphotoniso03",   &pfPhotonIso03oet );

  photonMVA_EE_2012->AddVariable("myphoton_sieie",   &sigietaieta );
  photonMVA_EE_2012->AddVariable("myphoton_sieip",   &sigietaiphi );
  photonMVA_EE_2012->AddVariable("myphoton_etawidth",   &scetawidth );
  photonMVA_EE_2012->AddVariable("myphoton_phiwidth",   &scphiwidth );
  photonMVA_EE_2012->AddVariable("myphoton_r9",   &r9 );
  
  photonMVA_EE_2012->AddVariable("myphoton_s4ratio",   &s4Ratio );
  
  photonMVA_EE_2012->AddVariable("myphoton_SCeta",   &etasc );
  photonMVA_EE_2012->AddVariable("event_rho",   &rho );
  
  photonMVA_EE_2012->AddVariable("myphoton_ESEffSigmaRR",   &sigRR);

  // photonMVA_EB_2012->BookMVA( methodName_Id, weightFile_IdEB_2012); //RAZOR UNNCESSARY
  // photonMVA_EE_2012->BookMVA( methodName_Id, weightFile_IdEE_2012);
  
  }

float HggPhotonID::getIdMVA(VecbosPho* pho, int nVertex, float rhoFastJet, int selVtxIndex){
  if(selVtxIndex < 0 || selVtxIndex >= vertices.size()){
    cout << "WARNING: Selected Vertex Index out of range: " << selVtxIndex << "/" << vertices.size() <<endl;
    return -9999;
  }
  if(pho->index <0) return -9999;
  if(debugPhotonID) cout << "Filling Variables" << endl;
  this->fillVariables(pho,nVertex,rhoFastJet,selVtxIndex);
  
  if(debugPhotonID) cout << "Preselection" << endl;
  if(! this->getPreSelection(pho,nVertex,rhoFastJet,selVtxIndex) ) return -9999;

  if(debugPhotonID) cout << "Version: " << version << endl;
  if(version.compare("May2012")==0) return (pho->isBarrel() ? photonMVA_EB_2012->EvaluateMVA(methodName_Id) : photonMVA_EE_2012->EvaluateMVA(methodName_Id) );    
  else return (pho->isBarrel() ? photonMVA_EB_2011->EvaluateMVA(methodName_Id) : photonMVA_EE_2011->EvaluateMVA(methodName_Id) );    
}

bool HggPhotonID::getIdCiC(VecbosPho* pho, int nVertex, float rhoFastJet,int selVtxIndex){
  if(selVtxIndex < 0 || selVtxIndex >= vertices.size()){
    cout << "WARNING: Selected Vertex Index out of range: " << selVtxIndex << "/" << vertices.size() <<endl;
    return false;
  }
  this->fillVariables(pho,nVertex,rhoFastJet,selVtxIndex);
  if(! this->getPreSelection(pho,nVertex,rhoFastJet,selVtxIndex) ) return false;
  //2011 CiC Cuts
  const int nCats=4;
  float cut_r9[nCats]         = {0.94,0.36,0.94,0.32};
  float cut_hoe[nCats]        = {0.082, 0.062,0.065,0.048};
  float cut_sieie[nCats]      = {0.0106,0.0097,0.028,0.027};
  float cut_combIso[nCats]    = {3.8,2.2,1.77,1.29};
  float cut_combIsoBad[nCats] = {11.7,3.4,3.9,1.84};
  float cut_trackIso[nCats]   = {3.5,2.2,2.3,1.45};

  int cat = this->getCiCCat(pho);

  if(pho->SC.r9 < cut_r9[cat]) return false;
  if(pho->HoverE > cut_hoe[cat]) return false;
  if(pho->SC.sigmaIEtaIEta > cut_sieie[cat]) return false;
  if(trkisooet > cut_trackIso[cat]) return false;
  if(isosumoet > cut_combIso[cat]) return false;
  if(isosumoetbad > cut_combIsoBad[cat]) return false;

  return true;
}

bool HggPhotonID::getIdCiCPF(VecbosPho* pho, int nVertex, float rhoFastJet,int selVtxIndex){
  if(selVtxIndex < 0 || selVtxIndex >= vertices.size()){
    cout << "WARNING: Selected Vertex Index out of range: " << selVtxIndex << "/" << vertices.size() <<endl;
    return false;
  }
  this->fillVariables(pho,nVertex,rhoFastJet,selVtxIndex);
  if(! this->getPreSelection(pho,nVertex,rhoFastJet,selVtxIndex) ) return false;
  //2012 CiC Cuts:
  const int nCats=4;
  float cut_r9[nCats]       = {0.94,0.298,0.94,0.24};
  float cut_hoe[nCats]      = {0.124, 0.092,0.142,0.063};
  float cut_sieie[nCats]    = {0.0108,0.0102,0.028,0.028};
  float cut_pfiso[nCats]    = {3.8,2.5,3.1,2.2};
  float cut_isoGood[nCats]  = {6.,4.7,5.6,3.6};
  float cut_isoBad[nCats]   = {10.,6.5,5.6,4.4};

  int cat = this->getCiCCat(pho);
  if(debugPhotonID) std::cout << "CiC Cat: " << cat << "   Photon SC Eta: " << pho->SC.eta << "  r9: " << pho->SC.r9 <<std::endl;

  if(debugPhotonID) std::cout << "H/E: " << pho->HoverE << "\tsieie: " << pho->SC.sigmaIEtaIEta
			      << "\tPFCharged: " << pfChargedIsoGood03oet << std::endl
			      << "isosum: " << isosumoetPF << "\tisosumoetbadPF: " << isosumoetbadPF << std::endl;
  if(pho->SC.r9 < cut_r9[cat]) return false;
  if(pho->HoverE > cut_hoe[cat]) return false;
  if(pho->SC.sigmaIEtaIEta > cut_sieie[cat]) return false;
  if(pfChargedIsoGood03oet > cut_pfiso[cat]) return false;
  if(isosumoetPF > cut_isoGood[cat]) return false;
  if(isosumoetbadPF > cut_isoBad[cat]) return false;

  return true;
}

bool HggPhotonID::getEGLooseID(VecbosPho* pho, int nVertex, float rhoFastJet,int selVtxIndex)
{
  if(selVtxIndex < 0 || selVtxIndex >= vertices.size()){
    cout << "WARNING: Selected Vertex Index out of range: " << selVtxIndex << "/" << vertices.size() <<endl;
    return false;
o  }

  //fill the variables and apply the preselection
  this->fillVariables(pho,nVertex,rhoFastJet,selVtxIndex);
  if(! this->getPreSelection(pho,nVertex,rhoFastJet,selVtxIndex) ) return false;

  //cuts for the barrel
  float sietaieta  = .012;
  float charged_had_iso03 = 2.6;
  float neutral_iso03 =  3.5; 
  float lin_neutral_iso03 = .04;
  float photon_iso03 = 1.3; 
  float lin_photon_iso03 = .005;
  float hoe = .05;

  //switch values if the photon is in the endcaps
  if(fabs(pho->SC.eta) > 1.48) {
    sietaieta = .034;
    charged_had_iso03 = 2.3;
    neutral_iso03 = 2.9;
    photon_iso03 = 999999999.9;
  }

  //calculate the rho corrections for the isolation
  float chosen_EA = 0;

  const int nCats = 7;
  float  EA_eta[nCats+1] = {0, 1, 1.48, 2.0, 2.2, 2.3, 2.4, 10000};
  float  EA_pho[nCats] = {0.148, 0.130, 0.112, 0.216, 0.262, 0.260, 0.266};

  for(int ii = 0; ii <= nCats; ii++) {
    float eta = pho->SC.eta;
    if (eta > EA_eta[ii] && eta < EA_eta[ii+1]){
      chosen_EA = [ii];
      break;
    }
  }

  float rhoCorr = chosen_EA * rhoFastJet;
  
  if(pho->HoverE > hoe) return false;
  if(pho->SC.sigmaIEtaIEta > sietaieta) return false;
  if((pfChargedIsoGood03 - rhoCorr)  > charged_had_iso03 ) return false;
  if((pho->dr03NeutralHadronPFIso - rhoCorr) > (neutral_iso03 + lin_neutral_iso03 * eT)) 
    return false;
  if( (pho->dr03PhotonPFIso - rhoCorr) > (photon_iso03 + lin_photon_iso03 * eT)) 
    return false;

  return true;
}

bool HggPhotonID::getIdCiCPF_Fake(VecbosPho* pho, int nVertex, float rhoFastJet,int selVtxIndex){
  if(selVtxIndex < 0 || selVtxIndex >= vertices.size()){
    cout << "WARNING: Selected Vertex Index out of range: " << selVtxIndex << "/" << vertices.size() <<endl;
    return false;
  }
  this->fillVariables(pho,nVertex,rhoFastJet,selVtxIndex);
  if(! this->getPreSelection(pho,nVertex,rhoFastJet,selVtxIndex) ) return false;
  //2012 CiC Cuts:
  const int nCats=4;
  float cut_r9[nCats]       = {0.94,0.298,0.94,0.24};
  float cut_hoe[nCats]      = {0.124, 0.092,0.142,0.063};
  float cut_sieie[nCats]    = {0.0108,0.0102,0.028,0.028};
  float cut_pfiso[nCats]    = {3.8,2.5,3.1,2.2};
  float cut_isoGood[nCats]  = {6.,4.7,5.6,3.6};
  float cut_isoBad[nCats]   = {10.,6.5,5.6,4.4};

  int cat = this->getCiCCat(pho);
  if(debugPhotonID) std::cout << "CiC Cat: " << cat << "   Photon SC Eta: " << pho->SC.eta << "  r9: " << pho->SC.r9 <<std::endl;

  if(debugPhotonID) std::cout << "H/E: " << pho->HoverE << "\tsieie: " << pho->SC.sigmaIEtaIEta
			      << "\tPFCharged: " << pfChargedIsoGood03oet << std::endl
			      << "isosum: " << isosumoetPF << "\tisosumoetbadPF: " << isosumoetbadPF << std::endl;
  
  //FAKE ID: USUAL CUTS, EXCEPT (!SIEIE  | !(2 COMB ISOLATIONS))
  
  if(pho->SC.r9 < cut_r9[cat]) return false;
  if(pho->HoverE > cut_hoe[cat]) return false;

  bool pass_sietaieta = pho->SC.sigmaIEtaIEta < cut_sieie[cat];
  bool pass_charged = pfChargedIsoGood03oet < cut_pfiso[cat];
  bool pass_sumgood = isosumoetPF < cut_isoGood[cat];
  bool pass_sumbad = isosumoetbadPF < cut_isoBad[cat];

  //catches all two fail scenarios except all passing
  bool one_pass = (pass_charged ^ pass_sumgood) ^ pass_sumbad;
  //catch the all passing scenario
  bool all_pass = pass_charged && pass_sumgood && pass_sumbad;

  //put top end limits on the isolation inversion
  bool no_high_iso = isosumoetPF < 20 && isosumoetbadPF < 100 && pfChargedIsoGood03oet < 20;
  bool pass_iso = one_pass && !all_pass;

  bool is_fake = (!pass_iso || !pass_sietaieta) && no_high_iso; 

  return is_fake;    
}

int HggPhotonID::getCiCCat(VecbosPho* pho){
  return (pho->SC.r9<0.94)+2*(fabs(pho->SC.eta) > 1.48);
}

void HggPhotonID::fillIsoVariables(VecbosPho* pho, ReducedPhotonData* data,int nVertex, float rhoFastJet,int selVtxIndex){
  this->fillVariables(pho,nVertex,rhoFastJet,selVtxIndex);
  
  data->HoverE = hoe;
  data->sieie  = sigietaieta;
  data->dr03PFChargedIso = pfChargedIsoGood03oet;
  data->isosumGood = isosumoetPF;
  data->isosumBad  = isosumoetbadPF;

  data->dr03EcalIso  = pho->dr03EcalRecHitSumEtCone - 0.012*eT;
  data->dr04HcalIso  = pho->dr03HcalTowerSumEtCone - 0.005*eT;
  data->dr03TrackIso = pho->dr03TrackIso[selVtxIndex] - 0.002*eT;
  data->dr02PFChargedIso = pho->dr02ChargedHadronPFIso[selVtxIndex];
}

bool HggPhotonID::getPreSelection(VecbosPho* pho, int nVertex, float rhoFastJet,int selVtxIndex){
  this->fillVariables(pho,nVertex,rhoFastJet,selVtxIndex);
  if(version.compare("May2012")==0) return this->getPreSelectionMay2012(pho,nVertex,rhoFastJet,selVtxIndex);
  else if(version.compare("Razor2013")==0) return this->getPreSelectionMay2012(pho,nVertex,rhoFastJet,selVtxIndex);
  else return this->getPreSelection2011(pho,nVertex,rhoFastJet,selVtxIndex);
}

bool HggPhotonID::getPreSelectionMay2012(VecbosPho* pho, int nVertex, float rhoFastJet, int selVtxIndex){
  if(debugPhotonID) std::cout << "getPreSelectionMay2012" << std::endl;

  if(debugPhotonID) std::cout << "et: " << eT << "\tecalIso: " << pho->dr03EcalRecHitSumEtCone -0.012*eT
			      << "\thcalIso: " << pho->dr03HcalTowerSumEtCone - 0.005*eT 
			      << "\ttrkIso:  " << pho->dr03TrackIso[selVtxIndex] - 0.002*eT
			      << "\tpfCharged: " << pho->dr02ChargedHadronPFIso[selVtxIndex] 
			      << std::endl;
  double ECALIso = (pho->dr03EcalRecHitSumEtCone - 0.012*eT);
  if(!doECALIso) ECALIso=0; // don't do th ECAL isolation

  if(pho->SC.r9 < 0.9){
    if( ECALIso > 4
        || (pho->dr03HcalTowerSumEtCone - 0.005*eT) > 4
        || (pho->dr03TrackIso[selVtxIndex] - 0.002*eT) > 4
	|| pho->dr02ChargedHadronPFIso[selVtxIndex] > 4) return false;
    if( (pho->isBarrel() && (pho->HoverE > 0.075 || pho->SC.sigmaIEtaIEta > 0.014) )
        || (!pho->isBarrel() && (pho->HoverE > 0.075 || pho->SC.sigmaIEtaIEta > 0.034) ) ) return false;
  }else{ //(SC->r9 > 0.9
    if( ECALIso > 50 
        || (pho->dr03HcalTowerSumEtCone - 0.005*eT) > 50
        || (pho->dr03TrackIso[selVtxIndex] - 0.002*eT) > 50 
	|| pho->dr02ChargedHadronPFIso[selVtxIndex] > 4) return false;
    if( (pho->isBarrel() && (pho->HoverE > 0.082 || pho->SC.sigmaIEtaIEta > 0.014) )
        || (!pho->isBarrel() && (pho->HoverE > 0.075 || pho->SC.sigmaIEtaIEta > 0.034) ) ) return false;
  }


  return true;
}

//deprecated
bool HggPhotonID::getPreSelectionRazor2013(VecbosPho* pho, int nVertex, float rhoFastJet, int selVtxIndex){
  if(debugPhotonID) std::cout << "getPreSelectionRazor2013" << std::endl;

  if(debugPhotonID) std::cout << "et: " << eT << "\tecalIso: " << pho->dr03EcalRecHitSumEtCone -0.012*eT
			      << "\thcalIso: " << pho->dr03HcalTowerSumEtCone - 0.005*eT 
			      << "\ttrkIso:  " << pho->dr03TrackIso[selVtxIndex] - 0.002*eT
			      << "\tpfCharged: " << pho->dr02ChargedHadronPFIso[selVtxIndex] 
			      << std::endl;
  double ECALIso = (pho->dr03EcalRecHitSumEtCone - 0.012*eT);
  if(!doECALIso) ECALIso=0; // don't do th ECAL isolation

  if(pho->SC.r9 < 0.9){
    if( ECALIso > 4
        || (pho->dr03HcalTowerSumEtCone - 0.005*eT) > 4
        || (pho->dr03TrackIso[selVtxIndex] - 0.002*eT) > 4
	|| pho->dr02ChargedHadronPFIso[selVtxIndex] > 4) return false;
    if( (pho->isBarrel() && (pho->HoverE > 0.1 || pho->SC.sigmaIEtaIEta > 0.014) )
        || (!pho->isBarrel() && (pho->HoverE > 0.1 || pho->SC.sigmaIEtaIEta > 0.034) ) ) return false;
  }else{ //(SC->r9 > 0.9
    if( ECALIso > 50 
        || (pho->dr03HcalTowerSumEtCone - 0.005*eT) > 50
        || (pho->dr03TrackIso[selVtxIndex] - 0.002*eT) > 50 
	|| pho->dr02ChargedHadronPFIso[selVtxIndex] > 4) return false;
    if( (pho->isBarrel() && (pho->HoverE > 0.1 || pho->SC.sigmaIEtaIEta > 0.014) )
        || (!pho->isBarrel() && (pho->HoverE > 0.1 || pho->SC.sigmaIEtaIEta > 0.034) ) ) return false;
  }


  return true;
}

bool HggPhotonID::getPreSelection2011(VecbosPho* pho, int nVertex, float rhoFastJet, int selVtxIndex){
  if(debugPhotonID) std::cout << "getPreSelection2011" << std::endl;

  if(pho->SC.r9 < 0.9){
    if( (pho->dr03EcalRecHitSumEtCone - 0.012*eT) > 4
        || (pho->dr04HcalTowerSumEtCone - 0.005*eT) > 4
        || (pho->dr03TrackIso[selVtxIndex] - 0.002*eT) > 4
        || pho->dr03EcalRecHitSumEtCone > 3
        || pho->dr04HcalTowerSumEtCone >  3
        || isosumoet > 2.8
        || pho->dr03TrkSumPtHollowCone > 4) return false;
    if( (pho->isBarrel() && (pho->HoverE > 0.075 || pho->SC.sigmaIEtaIEta > 0.014) )
        || (!pho->isBarrel() && (pho->HoverE > 0.075 || pho->SC.sigmaIEtaIEta > 0.034) ) ) return false;
  }else{ //(SC->r9 > 0.9
    if( (pho->dr03EcalRecHitSumEtCone - 0.012*eT) > 50
        || (pho->dr04HcalTowerSumEtCone - 0.005*eT) > 50
        || (pho->dr03TrackIso[selVtxIndex] - 0.002*eT) > 50
        || pho->dr03EcalRecHitSumEtCone > 3
        || pho->dr04HcalTowerSumEtCone >  3
        || isosumoet > 2.8
        || pho->dr03TrkSumPtHollowCone > 4) return false;
    if( (pho->isBarrel() && (pho->HoverE > 0.075 || pho->SC.sigmaIEtaIEta > 0.014) )
        || (!pho->isBarrel() && (pho->HoverE > 0.075 || pho->SC.sigmaIEtaIEta > 0.034) ) ) return false;
  }
  return true;
}
