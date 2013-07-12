#include <HggSelector.hh>
#include "ReadConfig.hh"
<<<<<<< HEAD
=======
#include <iostream>
#include <string>
#include <vector>
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a

//includes for the TMVA ID
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TRandom3.h"
#include "HggPhysUtils.cc"
using namespace std;
using namespace TMVA;

#define debugSelector 0

<<<<<<< HEAD
=======
//events
struct HggSelector::EventIndex {
  int RunNumber;
  Long64_t EventNumber;
  
  EventIndex() : RunNumber(0), EventNumber(0) {}
  
  bool operator <(const EventIndex &other) const
  {
    if(RunNumber < other.RunNumber)  return true;
    if(RunNumber > other.RunNumber)  return false;
    if(EventNumber < other.EventNumber) return true;
    if(EventNumber > other.EventNumber) return false;
    
    return false;
  }
};

>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
HggSelector::HggSelector():
  fChain(0),
  valid(false),
  doElectronVeto(true),
  doMuMuGamma(false),
  forceVtxZero(false),
  isData_(true),
  nSigma(3),
  Photons_(0),
  ggVerticesPhotonIndices(0),
  ggVerticesVertexIndex(0),
  ggVerticesPerEvtMVA(0),
  Muons_(0),
  Jets_(0),
  GenHiggs(0),
  GenPhotons(0)
{
}

HggSelector::HggSelector(vector<string> fNames, string treeName,string outFName):
  fChain(0),
  valid(false),
  doElectronVeto(true),
  doMuMuGamma(false),
  forceVtxZero(false),
  isData_(true),
  nSigma(3),
  Photons_(0),
  ggVerticesPhotonIndices(0),
  ggVerticesVertexIndex(0),
  ggVerticesPerEvtMVA(0),
  Muons_(0),
  Jets_(0),
  GenHiggs(0),
  GenPhotons(0)
{
  this->loadChain(fNames,treeName);
  outputFile = outFName;
}

HggSelector::~HggSelector(){
  delete fChain;
}


void HggSelector::loadChain(vector<string> fNames,string treeName){
  fChain = new TChain(treeName.c_str());
  vector<string>::const_iterator name;
  for(name = fNames.begin();name!=fNames.end();name++){
    fChain->AddFile(name->c_str());
  }
  valid = true;
}

int HggSelector::init(){
  if(!valid) return -1;

  mPair_ = -1;
  mPairNoCorr_ = -1;
  diPhoMVA_ = -999.;

  triggerDec = new int[triggers.size()];

<<<<<<< HEAD
=======
  //initialize the bad event list
  //  InitEventFlag("/home/jhardenbrook/2013/RAZOR_DIPHOTON/HggApp_Razor/AllBadABCDNEWTAUID.txt");
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
  
  this->setBranchAddresses();
  this->setupOutputTree();
  ReadConfig cfg;
  int retcode = cfg.read(configFile);
  if(retcode != 0){
    cout << "Error reading configuration file!" <<std::endl
	 << "Error Code: " << retcode << std::endl;
    valid = false;
    return -1;
  }
  weightFile_diPho = cfg.getParameter("weightFile_diPho");
  methodName_diPho  = cfg.getParameter("methodName_diPho");
  triggers = cfg.getTokens("Triggers",",");
  doRegression = cfg.getParameter("redoRegression").compare("yes")==0;
  doScale      = cfg.getParameter("redoScale").compare("yes")==0;
  doSmear      = cfg.getParameter("redoSmear").compare("yes")==0;

  if(doRegression){
    //corrector = new HggEGEnergyCorrector(this,
  }
  if(doScale){
    std::string scaleCFG = cfg.getParameter("EnergyScaleCFG");
    scale = new HggEnergyScale(scaleCFG);
  }
  if(doSmear){
    std::string smearCFG = cfg.getParameter("EnergySmearCFG");
    smear = new HggEnergyScale(smearCFG);    
    applyScaleSmear = atoi(cfg.getParameter("ScaleSmear").c_str());
    std::cout << "Doing smearing.  Config: " << smearCFG << "  ScaleSmear: " << applyScaleSmear <<std::endl;
  }

  cout << "Parameters: " << endl
       << weightFile_diPho << endl
       << methodName_diPho << endl;

  string MassResConfig = cfg.getParameter("MassResolutionConfig");
  massRes = new HggMassResolution();
  massRes->setConfigFile(MassResConfig);
  massRes->init();

  //get kinematic cuts
  string leadPhoEtMinS  =    cfg.getParameter("leadPhoEtMin");
  string subleadPhoEtMinS  = cfg.getParameter("subleadPhoEtMin");
  doPtOverM = !(cfg.getParameter("doPtOverM").compare("no")==0);
  string PtOverMLeadS =    cfg.getParameter("PtOverMLead");
  string PtOverMSubLeadS  = cfg.getParameter("PtOverMSubLead");

  leadPhoEtMin = 33.;
  subleadPhoEtMin = 25.;
  PtOverMLead = 40./120.;
  PtOverMSubLead = 30./120.;
  if(leadPhoEtMinS.compare("") != 0) leadPhoEtMin = atof(leadPhoEtMinS.c_str());
  if(subleadPhoEtMinS.compare("") != 0) subleadPhoEtMin = atof(subleadPhoEtMinS.c_str());
  if(PtOverMLeadS.compare("") != 0) PtOverMLead = atof(PtOverMLeadS.c_str());
  if(PtOverMSubLeadS.compare("") != 0) PtOverMSubLead = atof(PtOverMSubLeadS.c_str());

  cout << "Kinematic Selections: " << endl 
       << "pt Lead:   " << leadPhoEtMin << "   pt SubLead:   " << subleadPhoEtMin << endl
       << "Do Pt Over M: " << doPtOverM << endl
       << "pt/m Lead: " << PtOverMLead <<  "   pt/m SubLead: " << PtOverMSubLead << endl;

  PhotonID = new HggPhotonID();
  PhotonID->setConfig(configFile);
  PhotonID->Init();

  //setup MVA inputs histograms
  MVAInputs["massRes"] =  new TH1F("massRes","",200,0.,0.2);
  MVAInputs["massResWrongVtx"] = new TH1F("massResWrongVtx","",200,0.,0.2);
  MVAInputs["vtxProb"] = new TH1F("vtxProb","",100,0,1);
  MVAInputs["p1_EtByM"] = new TH1F("p1_EtByM","",200,0,1);
  MVAInputs["p2_EtByM"] = new TH1F("p2_EtByM","",200,0,1);
  MVAInputs["p1_Eta"] = new TH1F("p1_Eta","",120,-3,3);
  MVAInputs["p2_Eta"] = new TH1F("p2_Eta","",120,-3,3);
  MVAInputs["CosDPhi"] = new TH1F("cosDPhi","",200,-1,1);
  MVAInputs["p1_idMVA"] = new TH1F("p1_idMVA","",200,-1,1);
  MVAInputs["p2_idMVA"] = new TH1F("p2_idMVA","",200,-1,1);

  MjjDists["BeforeVBF"] = new TH1F("MjjBeforeVBF","",200,0,2000);
  MjjDists["dEta"] = new TH1F("dEta","",120,0,12);
  MjjDists["AfterDEta"] = new TH1F("MjjAfterDEta","",200,0,2000);
  MjjDists["Z"] = new TH1F("Z","",200,-20,20);
  MjjDists["AfterZ"] = new TH1F("MjjAfterZ","",200,0,2000);
  MjjDists["dPhi_jj_gg"] = new TH1F("dPhi_jj_gg","",32,0,3.2);
  MjjDists["Final"] = new TH1F("Final","",200,0,2000);


  this->setupTMVA();
  return 0;
}

void HggSelector::Loop(){
  if(!valid) return;

  if(this->init()!=0){
    std::cout << "ERROR INITIALIZING ... ABORTING!" <<std::endl;
    return;
  }
  this->setDefaults();
  cout << "Getting Entries ... "  << endl;
  Long64_t nEntries = fChain->GetEntries();
  Long64_t jentry=-1;
  int index1=-1,index2=-1;
  int index1PFCiC=-1,index2PFCiC=-1;
  int index1CiC=-1,index2CiC=-1;
<<<<<<< HEAD
=======

>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
  while(fChain->GetEntry(++jentry)){
    if(jentry%500==0) cout << ">> Processing Entry " << jentry << "/" << nEntries << endl;

    nPU_ = inPU;
    nVtxOut = nVtx;
    this->clear();
    if(!isData_) this->fillGenInfo();
    if(debugSelector) cout << "Setting Photon ID Vertices " << endl;
    PhotonID->setVertices(nVtx,vtxX,vtxY,vtxZ);

    if(debugSelector) cout << "requiring triggers ... " << endl;
    trigger_ = this->requireTrigger();

    // give the photon ID maker the vertex collection for this event

    if(debugSelector) cout << "Starting Photon Loop:  " << nPho_ << endl;
    for(int iPho=0; iPho<nPho_;iPho++){ //redo the scaling and smearing if needed
      VecbosPho *pho = &(Photons_->at(iPho));
      if(doScale && isData_){
	std::pair<float,float> dE = scale->getDEoE(*pho,runNumber);
	pho->dEoE    = dE.first;
        pho->dEoEErr = 0;
        pho->scaledEnergy = pho->correctedEnergy*(pho->dEoE);
        pho->scaledEnergyError = pho->correctedEnergyError*((pho->dEoE+pho->dEoEErr));
      }
      if(!isData_){
	pho->scaledEnergy = pho->correctedEnergy;
        pho->scaledEnergyError = pho->correctedEnergyError;
	if(doSmear){
	  std::pair<float,float> dE = smear->getDEoE(*pho,applyScaleSmear);
	  pho->dEoE    = dE.first;
	  pho->dEoEErr = dE.second;
	}
	smearPhoton(pho,0);	
      }
    }
    if(doMuMuGamma) this->fillMuMuGamma();
    

    if(debugSelector) cout << "done" << endl;
    if(debugSelector) cout << "# Photons: " << nPho_ << endl;

    if(debugSelector) cout << "Starting Photon Loop" << endl;

    float mvaOut[3];
    std::pair<int,int> indices;
    
    if(nSigma>0 && !isData_){ //do the energy scale and energy resolution systematics
      for(int iSmear=-nSigma;iSmear<=nSigma;iSmear++){
<<<<<<< HEAD
	//Do the Smearing for the MVA analysis
	indices = getBestPair(mvaOut,iSmear,0);
	diPhoMVASmear.push_back(mvaOut[0]);
	pho1MVASmear.push_back(mvaOut[1]);
	pho2MVASmear.push_back(mvaOut[2]);

	if(indices.first == -1 || indices.second==-1){
	  mPairSmear.push_back(-1);
	}else{
	  mPairSmear.push_back(getMPair(indices.first,indices.second));
	}
	//do the smearing for the PFCiC Analysis
	indices = getBestPairCiC(iSmear,0,true);
	if(indices.first == -1 || indices.second==-1){
	  mPairSmearPFCiC.push_back(-1);
	}else{
	  mPairSmearPFCiC.push_back(getMPair(indices.first,indices.second));
	}
	//do the smearing for the CiC Analysis
	indices = getBestPairCiC(iSmear,0,false);
	if(indices.first == -1 || indices.second==-1){
	  mPairSmearCiC.push_back(-1);
	}else{
	  mPairSmearCiC.push_back(getMPair(indices.first,indices.second));
	}
      } // Done with smearing

      for(int iScale=-nSigma;iScale<=nSigma;iScale++){ //do the scaling systematic
	indices = getBestPair(mvaOut,-999,iScale);
	diPhoMVAScale.push_back(mvaOut[0]);
	pho1MVAScale.push_back(mvaOut[1]);
	pho2MVAScale.push_back(mvaOut[2]);

	if(indices.first == -1 || indices.second==-1){
	  mPairScale.push_back(-1);
	}else{
	  mPairScale.push_back(getMPair(indices.first,indices.second));
	}
		//do the smearing for the PFCiC Analysis
	indices = getBestPairCiC(-999,iScale,true);
	if(indices.first == -1 || indices.second==-1){
	  mPairScalePFCiC.push_back(-1);
	}else{
	  mPairScalePFCiC.push_back(getMPair(indices.first,indices.second));
	}

		//do the smearing for the CiC Analysis
	indices = getBestPairCiC(-999,iScale,false);
	if(indices.first == -1 || indices.second==-1){
	  mPairScaleCiC.push_back(-1);
	}else{
	  mPairScaleCiC.push_back(getMPair(indices.first,indices.second));
	}
=======
        //Do the Smearing for the MVA analysis
        indices = getBestPair(mvaOut,iSmear,0);
        diPhoMVASmear.push_back(mvaOut[0]);
        pho1MVASmear.push_back(mvaOut[1]);
        pho2MVASmear.push_back(mvaOut[2]);

        if(indices.first == -1 || indices.second==-1){
          mPairSmear.push_back(-1);
        }else{
          mPairSmear.push_back(getMPair(indices.first,indices.second));
        }
        //do the smearing for the PFCiC Analysis
        indices = getBestPairCiC(iSmear,0,true);
        if(indices.first == -1 || indices.second==-1){
          mPairSmearPFCiC.push_back(-1);
        }else{
          mPairSmearPFCiC.push_back(getMPair(indices.first,indices.second));
        }
        //do the smearing for the CiC Analysis
        indices = getBestPairCiC(iSmear,0,false);
        if(indices.first == -1 || indices.second==-1){
          mPairSmearCiC.push_back(-1);
        }else{
          mPairSmearCiC.push_back(getMPair(indices.first,indices.second));
        }
      } // Done with smearing
      
      for(int iScale=-nSigma;iScale<=nSigma;iScale++){ //do the scaling systematic
        indices = getBestPair(mvaOut,-999,iScale);
        diPhoMVAScale.push_back(mvaOut[0]);
        pho1MVAScale.push_back(mvaOut[1]);
        pho2MVAScale.push_back(mvaOut[2]);
        
        if(indices.first == -1 || indices.second==-1){
          mPairScale.push_back(-1);
        }else{
          mPairScale.push_back(getMPair(indices.first,indices.second));
        }
		//do the smearing for the PFCiC Analysis
        indices = getBestPairCiC(-999,iScale,true);
        if(indices.first == -1 || indices.second==-1){
          mPairScalePFCiC.push_back(-1);
        }else{
          mPairScalePFCiC.push_back(getMPair(indices.first,indices.second));
        }
        
		//do the smearing for the CiC Analysis
        indices = getBestPairCiC(-999,iScale,false);
        if(indices.first == -1 || indices.second==-1){
          mPairScaleCiC.push_back(-1);
        }else{
          mPairScaleCiC.push_back(getMPair(indices.first,indices.second));
        }
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
      }
    }
    
    indices = getBestPair(&diPhoMVA_,0,0); // no scaling and default smearing
    index1 = indices.first;
    index2 = indices.second;
    indices = make_pair(-1,-1);
    //Do the PFCiC selection as well!
    indices = getBestPairCiC(0,0,true); // no scaling and default smearing
    index1PFCiC = indices.first;
    index2PFCiC = indices.second;
    indices = make_pair(-1,-1);

    //Do the CiC selection as well!
    indices = getBestPairCiC(0,0,false); // no scaling and default smearing
    index1CiC = indices.first;
    index2CiC = indices.second;
    indices = make_pair(-1,-1);

    if(debugSelector) cout << "LOOP DONE" << endl;	  
    

    if(debugSelector) cout << "indices: " << index1 << "  " << index2 << endl;
    if(debugSelector) cout << "indicesPFCiC: " << index1PFCiC << "  " << index2PFCiC << endl;
    if(debugSelector) cout << "indicesCiC: " << index1CiC << "  " << index2CiC << endl;

<<<<<<< HEAD
    if(index1 > -1 && index2 > -1){
=======
    //    if(index1 > -1 && index2 > -1){
    if(false){ /////////DO NOT FILL MVA  VARIABLES FOR RAZOR ANALYSIS////////
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
      //fill MVA variables
      int selectedVertex = getVertexIndex(index1,index2);
      if(debugSelector) cout << "Final Selection MVA: " << selectedVertex << endl;
      
      TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
      pho1_ = Photons_->at(index1);
      pho2_ = Photons_->at(index2);
      OutPhotons_.push_back(getReducedData(&pho1_,vtxPos,selectedVertex));
      OutPhotons_.push_back(getReducedData(&pho2_,vtxPos,selectedVertex));

      if(debugSelector) cout << "Photon Indices: " << pho1_.index << "  " << pho2_.index << endl;

      mPair_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
      mPairNoCorr_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
      mPairRes_ = massRes->getMassResolutionEonly(&pho1_,&pho2_,vtxPos);
      mPairResWrongVtx_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
      diPhoVtx_ = selectedVertex;
      diPhoVtxX_ = vtxX[selectedVertex];
      diPhoVtxY_ = vtxY[selectedVertex];
      diPhoVtxZ_ = vtxZ[selectedVertex];
      vtxProb_ = this->getVertexProb(index1,index2);

      float jpt[2];
      Mjj_  = this->getVBFMjj(&pho1_,&pho2_,vtxPos,jpt);
      ptJet1_ = jpt[0];
      ptJet2_ = jpt[1];

      TLorentzVector p1 = pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy);
      TLorentzVector p2 = pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy);
      if(p1.Pt() < p2.Pt()){
	TLorentzVector tmp = p1;
	p1=p2; p2=tmp;
      }
      TLorentzVector gg = p1+p2;
      TVector3 boost = -1*gg.BoostVector();
      p1.Boost(boost);
      cosThetaLead = p1.Vect().Dot(gg.Vect())/p1.Vect().Mag()/gg.Vect().Mag();

      cat_ = getCategory();


      float mva1 = PhotonID->getIdMVA(&pho1_,nVtx,rho,selectedVertex);
      float mva2 = PhotonID->getIdMVA(&pho2_,nVtx,rho,selectedVertex);

      for(int i1=0;i1<3;i1++){ for(int i2=0;i2<3;i2++){
	  float offset1 = 0.01*(i1-1);
	  float offset2 = 0.01*(i2-1);
	  diPhoMVAShift_[3*i1+i2] =  getDiPhoMVA(index1,index2,mva1+offset1,mva2+offset2,false);
	} }



    }else{
      mPair_=-1;      
      cat_=-1;
    }
<<<<<<< HEAD

=======
    
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
    if(index1PFCiC > -1 && index2PFCiC > -1){
      //fill PFCiC variables
      int selectedVertex = getVertexIndex(index1PFCiC,index2PFCiC);
      if(debugSelector) cout << "Final Selection PFCiC: " << selectedVertex << endl;
      
      TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
      pho1_ = Photons_->at(index1PFCiC);
      pho2_ = Photons_->at(index2PFCiC);
      OutPhotonsPFCiC_.push_back(getReducedData(&pho1_,vtxPos,selectedVertex));
      OutPhotonsPFCiC_.push_back(getReducedData(&pho2_,vtxPos,selectedVertex));

      mPairPFCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
      mPairNoCorrPFCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
      mPairResPFCiC_ = massRes->getMassResolutionEonly(&pho1_,&pho2_,vtxPos);
      mPairResWrongVtxPFCiC_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
      diPhoVtxPFCiC_ = selectedVertex;
      diPhoVtxXPFCiC_ = vtxX[selectedVertex];
      diPhoVtxYPFCiC_ = vtxY[selectedVertex];
      diPhoVtxZPFCiC_ = vtxZ[selectedVertex];
      vtxProbPFCiC_ = this->getVertexProb(index1PFCiC,index2PFCiC);

      float jpt[2];
      MjjPFCiC_  = this->getVBFMjj(&pho1_,&pho2_,vtxPos,jpt);
      ptJet1PFCiC_ = jpt[0];
      ptJet2PFCiC_ = jpt[1];

      TLorentzVector p1 = pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy);
      TLorentzVector p2 = pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy);
      if(p1.Pt() < p2.Pt()){
<<<<<<< HEAD
	TLorentzVector tmp = p1;
	p1=p2; p2=tmp;
=======
        TLorentzVector tmp = p1;
        p1=p2; p2=tmp;
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
      }
      TLorentzVector gg = p1+p2;
      TVector3 boost = -1*gg.BoostVector();
      p1.Boost(boost);
      cosThetaLeadPFCiC = p1.Vect().Dot(gg.Vect())/p1.Vect().Mag()/gg.Vect().Mag();
<<<<<<< HEAD

      catPFCiC_ = getCategoryPFCiC();
    }else{
      mPairPFCiC_=-1;
      catPFCiC_ = -1;
    }

    if(index1CiC > -1 && index2CiC > -1){
=======
      
      catPFCiC_ = getCategoryPFCiC();
      
      //////////////CALCULATE RAZOR VARIABLES/////////////
      //////////////ONLY FOR PFCIC ANALYSIS//////////////      

      //get the jet list
      vector<TLorentzVector> jetlist = GetJetList(p1,p2);
      
      nJets = jetlist.size();

      //cut on the number of jets in the event
      float min_jet_cut = 1;


      //must be barrel photons (Use the supercluster)
      bothBarrel = fabs(pho1_.SC.eta) < 1.48 && fabs(pho2_.SC.eta) < 1.48;      
      //check the bad event list (needs to be produced)
      //      badEventList = isFlagged();
      bool bothEndcaps = fabs(pho1_.SC.eta) < 2.6 && fabs(pho2_.SC.eta) < 2.6;

      //must pass met filters! and not be in the bad event list
      bool calcRazor = (jetlist.size() >= min_jet_cut) && PassMETFilters() && bothEndcaps && !badEventList;

      if(calcRazor) {		
	//combine the photons and jets into hemispheres      
	vector<TLorentzVector> tmpJet = CombineJets_R_no_seed(jetlist, p1, p2);
	vector<TLorentzVector> tmpJet_SS = CombineJets_R_SSorOS(jetlist, p1, p2, true);
	vector<TLorentzVector> tmpJet_OS = CombineJets_R_SSorOS(jetlist, p1, p2, false);

	if((tmpJet.size() != 2) || (tmpJet.size() != 2) || (tmpJet.size() != 2)) {
	  cout << "HEMISPHERE MAKER DOES NOT RETURN 2 HEMISPHERES" << endl;
	}
	
        TLorentzVector PFHem1 = tmpJet[0];
        TLorentzVector PFHem2 = tmpJet[1];

        TLorentzVector PFHem1_SS = tmpJet_SS[0];
        TLorentzVector PFHem2_SS = tmpJet_SS[1];

        TLorentzVector PFHem1_OS = tmpJet_OS[0];
        TLorentzVector PFHem2_OS = tmpJet_OS[1];
        
        //calculate the variables, no seed, SS, and OS
        double MT = CalcMTR(PFHem1, PFHem2, pfMet);
        double MT_SS = CalcMTR(PFHem1_SS, PFHem2_SS, pfMet);
        double MT_OS = CalcMTR(PFHem1_OS, PFHem2_OS, pfMet);

        double variable = -999999.;
        double Rvariable = -999999.;
        double variable_SS = -999999.;
        double Rvariable_SS = -999999.;
        double variable_OS = -999999.;
        double Rvariable_OS = -999999.;

        variable = CalcGammaMRstar(PFHem1, PFHem2);
        variable_SS = CalcGammaMRstar(PFHem1_SS, PFHem2_SS);
        variable_OS = CalcGammaMRstar(PFHem1_OS, PFHem2_OS);

        if(variable > 0) Rvariable = MT/variable;
        if(variable_SS > 0) Rvariable_SS = MT_SS/variable_SS;
        if(variable_OS > 0) Rvariable_OS = MT_OS/variable_OS;
        
        //razor variables
        //no seeding
        PFMR = variable;
        PFR = Rvariable;
        //same side photons
        PFMR_SS = variable_SS;
        PFR_SS = Rvariable_SS;
        //opposite side
        PFMR_OS = variable_OS;
        PFR_OS = Rvariable_OS;

        nBtags = 0;
        
        //no seeding hemispheres
        ptHem1 = PFHem1.Pt();
        etaHem1 = PFHem1.Eta();
        phiHem1 =  PFHem1.Phi();
        
        ptHem2 = PFHem2.Pt();
        etaHem2 = PFHem2.Eta();
        phiHem2 = PFHem2.Phi();      

        mHem1 = PFHem1.E();
        mHem2 = PFHem2.E();

        //same side hemispheres
        ptHem1_SS = PFHem1_SS.Pt();
        etaHem1_SS = PFHem1_SS.Eta();
        phiHem1_SS =  PFHem1_SS.Phi();
        
        ptHem2_SS = PFHem2_SS.Pt();
        etaHem2_SS = PFHem2_SS.Eta();
        phiHem2_SS = PFHem2_SS.Phi();      

        mHem1_SS = PFHem1_SS.E();
        mHem2_SS = PFHem2_SS.E();

        //opposite side hemispheres
        ptHem1_OS = PFHem1_OS.Pt();
        etaHem1_OS = PFHem1_OS.Eta();
        phiHem1_OS =  PFHem1_OS.Phi();
        
        ptHem2_OS = PFHem2_OS.Pt();
        etaHem2_OS = PFHem2_OS.Eta();
        phiHem2_OS = PFHem2_OS.Phi();      

        mHem1_OS = PFHem1_OS.E();
        mHem2_OS = PFHem2_OS.E();

      }
      else{ // calcRazor == false
        FillRazorVarsWith(0);
      }
    }
    else{
      mPairPFCiC_=-1;
      catPFCiC_ = -1;
      //To reduce file size, only fill the tree if we reconstruct the razor
      continue;
      FillRazorVarsWith(-99);
    }

    //if(index1CiC > -1 && index2CiC > -1){
    if(false){ /////////DO NOT FILL CiC VARIABLES FOR RAZOR ANALYSIS////////
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
      //fill CiC variables
      int selectedVertex = getVertexIndex(index1CiC,index2CiC);
      if(debugSelector) cout << "Final Selection CiC: " << selectedVertex << endl;
      
      TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
      pho1_ = Photons_->at(index1CiC);
      pho2_ = Photons_->at(index2CiC);
      OutPhotonsCiC_.push_back(getReducedData(&pho1_,vtxPos,selectedVertex));
      OutPhotonsCiC_.push_back(getReducedData(&pho2_,vtxPos,selectedVertex));
      if(debugSelector) cout << "Done Getting Photons" << endl;

      mPairCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy) + pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy)).M();
      mPairNoCorrCiC_ = (pho1_.p4FromVtx(vtxPos,pho1_.energy) + pho2_.p4FromVtx(vtxPos,pho2_.energy)).M();
      mPairResCiC_ = massRes->getMassResolutionEonly(&pho1_,&pho2_,vtxPos);
      mPairResWrongVtxCiC_ = massRes->getMassResolution(&pho1_,&pho2_,vtxPos,true);
      diPhoVtxCiC_ = selectedVertex;
      diPhoVtxXCiC_ = vtxX[selectedVertex];
      diPhoVtxYCiC_ = vtxY[selectedVertex];
      diPhoVtxZCiC_ = vtxZ[selectedVertex];
      vtxProbCiC_ = this->getVertexProb(index1CiC,index2CiC);

      float jpt[2];
      MjjCiC_  = this->getVBFMjj(&pho1_,&pho2_,vtxPos,jpt);
      ptJet1CiC_ = jpt[0];
      ptJet2CiC_ = jpt[1];

      TLorentzVector p1 = pho1_.p4FromVtx(vtxPos,pho1_.finalEnergy);
      TLorentzVector p2 = pho2_.p4FromVtx(vtxPos,pho2_.finalEnergy);
      if(p1.Pt() < p2.Pt()){
	TLorentzVector tmp = p1;
	p1=p2; p2=tmp;
      }
      TLorentzVector gg = p1+p2;
      TVector3 boost = -1*gg.BoostVector();
      p1.Boost(boost);
      cosThetaLeadCiC = p1.Vect().Dot(gg.Vect())/p1.Vect().Mag()/gg.Vect().Mag();
    }else{
      mPairCiC_=-1;
    }

    nOutPhotons_ = OutPhotons_.size();
    nOutPhotonsPFCiC_ = OutPhotonsPFCiC_.size();
    nOutPhotonsCiC_ = OutPhotonsCiC_.size();

    MET = pfMet;
    METPhi = pfMetPhi;
    lumiBlockOut = lumiBlock;
    runNumberOut = runNumber;
    evtNumberOut = evtNumber;
    evtWeight = pileupWeight;

    eeBadScFilterFlagOut         = eeBadScFilterFlag;
    hcalLaserEventFilterFlagOut  = hcalLaserEventFilterFlag;
    HBHENoiseFilterResultFlagOut = HBHENoiseFilterResultFlag;
    isNotDeadEcalClusterOut      = isNotDeadEcalCluster;
    trackerFailureFilterFlagOut  = trackerFailureFilterFlag;
    CSCHaloFilterFlagOut         = CSCHaloFilterFlag;
    drDeadOut                    = drDead; 
    drBoundaryOut                = drBoundary;
    ECALTPFilterFlagOut          = ECALTPFilterFlag;

<<<<<<< HEAD
=======
    
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
    outTree->Fill();
  }//while(fChain...

TFile *f = new TFile(outputFile.c_str(),"RECREATE");
  outTree->Write();
  if(doMuMuGamma) outTreeMuMuG->Write();
  std::map<std::string,TH1F*>::iterator it;
  for(it = MVAInputs.begin(); it!=MVAInputs.end(); it++) (*it).second->Write();
  for(it = MjjDists.begin();  it!=MjjDists.end();  it++) (*it).second->Write();
  if(PhotonID->getHists()){
    for(it = PhotonID->getHists()->begin();  it!=PhotonID->getHists()->end();  it++) (*it).second->Write();

  }

  f->Close();
}

bool HggSelector::preSelectPhotons(VecbosPho* pho1,VecbosPho* pho2,TVector3 vtxPos){
  //apply kinematic photon selection
  if(fabs(pho1->SC.eta) > 2.5 || fabs(pho2->SC.eta) > 2.5) return false;  //outside of tracker acceptance
  if(fabs(pho1->SC.eta) > 1.4442 && fabs(pho1->SC.eta) < 1.566) return false;
  if(fabs(pho2->SC.eta) > 1.4442 && fabs(pho2->SC.eta) < 1.566) return false; // veto gap photons

  if( pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt() < subleadPhoEtMin || pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt() < subleadPhoEtMin ) return false;
  if( pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt() < leadPhoEtMin && pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt() < leadPhoEtMin ) return false;

  if(doPtOverM){
    float M = (pho1->p4FromVtx(vtxPos,pho1->finalEnergy) + pho2->p4FromVtx(vtxPos,pho2->finalEnergy)).M();
    if( pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()/M < PtOverMSubLead || pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt()/M < PtOverMSubLead ) return false;
    if( pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()/M < PtOverMLead && pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt()/M < PtOverMLead ) return false;
    
  }

  return true;
}

float HggSelector::getMPair(int i1, int i2){
  int selectedVertex = getVertexIndex(i1,i2);
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);
  
  VecbosPho pho1 = Photons_->at(i1);
  VecbosPho pho2 = Photons_->at(i2);
  
  return (pho1.p4FromVtx(vtxPos,pho1.finalEnergy) + pho2.p4FromVtx(vtxPos,pho2.finalEnergy)).M();
}

std::pair<int,int> HggSelector::getBestPairCiC(int smearShift,int scaleShift,bool usePF=true){
  std::pair<int,int> indices(-1,-1);
  //int bestCat=-1;
  //float bestMass=-99;
  float highestPtSum=0; // is this the best way to select the pairs?
    for(int iPho1=0; iPho1<nPho_;iPho1++){
      if(photonMatchedElectron[iPho1] && doElectronVeto) continue;
      for(int iPho2=iPho1; iPho2<nPho_;iPho2++){
<<<<<<< HEAD
	if(iPho1==iPho2) continue;
	if(photonMatchedElectron[iPho2] && doElectronVeto) continue;
	if(debugSelector) cout << ">> " << iPho1 << "  " << iPho2 << endl;
	//scale/smear the energy of the photon
	VecbosPho* pho1 = &(Photons_->at(iPho1));
	VecbosPho* pho2 = &(Photons_->at(iPho2));
	int selVtxI = this->getVertexIndex(iPho1,iPho2);
	TVector3 vtxPos(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]);
	if(!this->preSelectPhotons(pho1,pho2,vtxPos)) continue;
	if(!isData_){
	  
	  //apply scale shift	  
	  pho1->finalEnergy = pho1->scaledEnergy + scaleShift*pho1->scaledEnergyError;
	  pho2->finalEnergy = pho2->scaledEnergy + scaleShift*pho2->scaledEnergyError;

	  smearPhoton(pho1,smearShift);
	  smearPhoton(pho2,smearShift);
	}
	bool CiC1,CiC2;
	if(usePF){
	  CiC1 = PhotonID->getIdCiCPF(pho1,nVtx,rho,selVtxI);
	  CiC2 = PhotonID->getIdCiCPF(pho2,nVtx,rho,selVtxI);
	}else{
	  CiC1 = PhotonID->getIdCiC(pho1,nVtx,rho,selVtxI);
	  CiC2 = PhotonID->getIdCiC(pho2,nVtx,rho,selVtxI);
	}
	if(!CiC1 || !CiC2) continue;
	float thisPtSum = pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()
	  + pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt();	
	if(thisPtSum > highestPtSum){
	  highestPtSum = thisPtSum;
	  indices.first = iPho1;
	  indices.second = iPho2;
	}
=======
        if(iPho1==iPho2) continue;
        if(photonMatchedElectron[iPho2] && doElectronVeto) continue;
        if(debugSelector) cout << ">> " << iPho1 << "  " << iPho2 << endl;
        //scale/smear the energy of the photon
        VecbosPho* pho1 = &(Photons_->at(iPho1));
        VecbosPho* pho2 = &(Photons_->at(iPho2));
        int selVtxI = this->getVertexIndex(iPho1,iPho2);
        TVector3 vtxPos(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]);
        if(!this->preSelectPhotons(pho1,pho2,vtxPos)) continue;
        if(!isData_){
	  
          //apply scale shift	  
          pho1->finalEnergy = pho1->scaledEnergy + scaleShift*pho1->scaledEnergyError;
          pho2->finalEnergy = pho2->scaledEnergy + scaleShift*pho2->scaledEnergyError;
          
          smearPhoton(pho1,smearShift);
          smearPhoton(pho2,smearShift);
        }
        bool CiC1,CiC2;
        if(usePF){
          CiC1 = PhotonID->getIdCiCPF(pho1,nVtx,rho,selVtxI);
          CiC2 = PhotonID->getIdCiCPF(pho2,nVtx,rho,selVtxI);
        }else{
          CiC1 = PhotonID->getIdCiC(pho1,nVtx,rho,selVtxI);
          CiC2 = PhotonID->getIdCiC(pho2,nVtx,rho,selVtxI);
        }
        if(!CiC1 || !CiC2) continue;
        float thisPtSum = pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()
          + pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt();	
        if(thisPtSum > highestPtSum){
          highestPtSum = thisPtSum;
          indices.first = iPho1;
          indices.second = iPho2;
        }
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
      }// for(iPho2...
    }// for(iPho1...
    return indices;
}

std::pair<int,int> HggSelector::getBestPair(float* mvaOut, int smearShift,int scaleShift){
  std::pair<int,int> indices(-1,-1);
  float diPhoMVAMax=-99;
  float maxSumPt=-1;
  *mvaOut = -999.;
  for(int iPho1=0; iPho1<nPho_;iPho1++){
    if(photonMatchedElectron[iPho1] && doElectronVeto) continue;
    for(int iPho2=iPho1; iPho2<nPho_;iPho2++){
      if(iPho1==iPho2) continue;
      if(photonMatchedElectron[iPho2] && doElectronVeto) continue;
      if(debugSelector) cout << ">> " << iPho1 << "  " << iPho2 << endl;
      //NON-PF block
      //scale/smear the energy of the photon
      VecbosPho* pho1 = &(Photons_->at(iPho1));
      VecbosPho* pho2 = &(Photons_->at(iPho2));
      int selVtxI = this->getVertexIndex(iPho1,iPho2);
      TVector3 vtxPos(vtxX[selVtxI],vtxY[selVtxI],vtxZ[selVtxI]);
      if(!this->preSelectPhotons(pho1,pho2,vtxPos)) continue;
      if(!isData_){
	
<<<<<<< HEAD
	//apply scale shift	  
	pho1->finalEnergy = pho1->scaledEnergy + scaleShift*pho1->scaledEnergyError;
	pho2->finalEnergy = pho2->scaledEnergy + scaleShift*pho2->scaledEnergyError;
	smearPhoton(pho1,smearShift);
	smearPhoton(pho2,smearShift);
=======
        //apply scale shift	  
        pho1->finalEnergy = pho1->scaledEnergy + scaleShift*pho1->scaledEnergyError;
        pho2->finalEnergy = pho2->scaledEnergy + scaleShift*pho2->scaledEnergyError;
        smearPhoton(pho1,smearShift);
        smearPhoton(pho2,smearShift);
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
      }
      if(debugSelector) cout << "Getting Photon ID:" << endl;
      float mva1 = PhotonID->getIdMVA(pho1,nVtx,rho,selVtxI);
      if(debugSelector) cout << "Getting Photon ID:" << endl;
      float mva2 = PhotonID->getIdMVA(pho2,nVtx,rho,selVtxI);
      if(debugSelector) cout << "Getting Double Photon MVA" << endl;
      float diPhoMVA =  getDiPhoMVA(iPho1,iPho2,mva1,mva2,false);
      if(debugSelector) cout << "\t\t" << mva1 << "  " << mva2 << "  " << diPhoMVA << endl;
      float thisPtSum = pho1->p4FromVtx(vtxPos,pho1->finalEnergy).Pt()
<<<<<<< HEAD
	+ pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt();	

      if(diPhoMVA>=-1 &&  thisPtSum > maxSumPt){
	maxSumPt = thisPtSum;
	indices.first = iPho1;
	indices.second = iPho2;
	diPhoMVAMax = diPhoMVA;
	*mvaOut = diPhoMVA;
=======
        + pho2->p4FromVtx(vtxPos,pho2->finalEnergy).Pt();	
      
      if(diPhoMVA>=-1 &&  thisPtSum > maxSumPt){
        maxSumPt = thisPtSum;
        indices.first = iPho1;
        indices.second = iPho2;
        diPhoMVAMax = diPhoMVA;
        *mvaOut = diPhoMVA;
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
      }
    }// for(iPho2...
  }// for(iPho1...
  
    return indices;
}

void HggSelector::smearPhoton(VecbosPho* pho,int smearShift){
  TRandom3 rng(0);
  float smear = pho->dEoE + smearShift*pho->dEoEErr;
  if(smear < 0) smear = 0;
  if(smearShift<-100) smear = 0;
  float rand=0;
  if(smear > 0) rand = rng.Gaus(0,smear);
  if(rand < -1) rand=-1;
  if(rand > 1E3) rand = 1E3;
  pho->finalEnergy = pho->scaledEnergy*(1+rand);
  pho->finalEnergyError = pho->scaledEnergyError;//*(1+rand);
}

void HggSelector::setDefaults(){
  mPair_=-1;
  mPairNoCorr_=-1;
  genHiggsPt = -1;
  nPU_ = inPU;
  nVtxOut = nVtx;
  evtWeight=1;
}
void HggSelector::clear(){
  mPair_=-1;
  mPairNoCorr_=-1;
  mPairPFCiC_=-1;
  mPairNoCorrPFCiC_=-1;
  mPairCiC_=-1;
  mPairNoCorrCiC_=-1;
  diPhoMVA_=-999;
  diPhoVtx_= -1;
  diPhoVtxX_=0;
  diPhoVtxY_=0;
  diPhoVtxZ_=0;

  Mjj_ = 0;
  MjjPFCiC_ = 0;
  MjjCiC_ = 0;

  if(debugSelector) std::cout << "Clearing Scale/Smear variables" << std::endl;
  mPairScale.clear();
  pho1MVAScale.clear();
  pho2MVAScale.clear();
  diPhoMVAScale.clear();
  mPairSmear.clear();
  pho1MVASmear.clear();
  pho2MVASmear.clear();
  diPhoMVASmear.clear();

  mPairScalePFCiC.clear();
  mPairSmearPFCiC.clear();

  mPairScaleCiC.clear();
  mPairSmearCiC.clear();

  if(debugSelector) std::cout << "Clearing Output Variables" << std::endl;
  OutPhotons_.clear();
  OutPhotonsPFCiC_.clear();
  OutPhotonsCiC_.clear();

  if(debugSelector) std::cout << "Clearing MMG Output Variables" << std::endl;
  if(doMuMuGamma){
    MMG_Mu1.clear();
    MMG_Mu2.clear();
    MMG_Pho.clear();
  }
}

void HggSelector::fillGenInfo(){
  if(debugSelector) std::cout << "Setting Gen Higgs Info" << std::endl;
  if(nGenHiggs > 0){
    genHiggsPt = GenHiggs->front().pt;
    genHiggsVx = GenHiggs->front().Vx;
    genHiggsVy = GenHiggs->front().Vy;
    genHiggsVz = GenHiggs->front().Vz;
  }else{
    genHiggsPt = 0;
    genHiggsVx = 0;
    genHiggsVy = 0;
    genHiggsVz = 0;
  }

  if(debugSelector) std::cout << "Clearing Output Variables" << std::endl;
  int selGenPho=0;
  GenCollection::const_iterator genPho;
  if(GenPhotons == 0){
    std::cout << "WARNING: GenPhotons Collection is EMPTY!" << std::endl;
    return;
  }
  
  for(int i=0;i<nGenPho;i++){
    VecbosGen* genPho = &(GenPhotons->at(i));
    if(genPho->idMother != 25) continue; // not a higgs
    if(selGenPho==0){
      etaGenPho1 = genPho->eta;
      phiGenPho1 = genPho->phi;
      ptGenPho1 =  genPho->pt;
      energyGenPho1 = genPho->energy;
    }
    else if(selGenPho==1){
      etaGenPho2 = genPho->eta;
      phiGenPho2 = genPho->phi;
      ptGenPho2 =  genPho->pt;
      energyGenPho2 = genPho->energy;
    }else{
      cout << "More than 2 photons from a higgs decay!!" << endl;
    }
    selGenPho++;
  }
  
}

ReducedPhotonData HggSelector::getReducedData(VecbosPho* pho,TVector3 selVtx,int selVtxI){
  if(debugSelector) cout << "Filling Reduced Data: " << flush; 
  if(debugSelector) cout << pho->index << "  " << selVtx.Z() << "  " <<selVtxI <<endl;

  ReducedPhotonData data;
  TLorentzVector p4 = pho->p4FromVtx(selVtx,pho->finalEnergy,false);
  TLorentzVector p4NoCorr = pho->p4FromVtx(selVtx,pho->energy,false);
  TLorentzVector p4Gen;
  if(pho->genMatch.index!=-1) {
    p4Gen.SetPtEtaPhiE(pho->genMatch.pt,pho->genMatch.eta,pho->genMatch.phi,pho->genMatch.energy);
    data.pt_Gen = p4Gen.Pt(); data.eta_Gen = p4Gen.Eta(); data.phi_Gen = p4Gen.Phi(); data.E_Gen = p4Gen.E();
  }
  else{
    p4Gen.SetPtEtaPhiE(0.,0.,0.,0.);
    data.pt_Gen = 0.; data.eta_Gen = 0.; data.phi_Gen = 0.; data.E_Gen = 0.;  
  }
  data.pt = p4.Pt(); data.eta = p4.Eta(); data.phi = p4.Phi(); data.E =p4.E(); data.EError = pho->finalEnergyError;
  data.EErrorSmeared = massRes->getResolution(pho);
  data.index = pho->index;
  data.etaSC = pho->SC.eta;
  data.r9 = pho->SC.r9;
  data.passPFCiC = PhotonID->getIdCiCPF(pho,nVtx,rho,selVtxI); 
  data.category = (data.r9 < 0.94)+2*(fabs(data.etaSC) > 1.48); 
  data.idMVA = PhotonID->getIdMVA(pho,nVtx,rho,selVtxI);  
  data.mother = pho->genMatch.idMother; 
  PhotonID->fillIsoVariables(pho,&data,nVtx,rho,selVtxI);
  if(debugSelector) cout << "DONE Filling Reduced Data" <<endl;
  return data;
}

void HggSelector::setupTMVA(){
  diPhotonMVA = new TMVA::Reader( "!Color:!Silent" );; 

  // diphoton                
  diPhotonMVA->AddVariable("masserrsmeared/mass",&smearedMassErrByMass); 
  diPhotonMVA->AddVariable("masserrsmearedwrongvtx/mass",&smearedMassErrByMassWrongVtx);         
  diPhotonMVA->AddVariable("vtxprob",&vtxprob);             
  diPhotonMVA->AddVariable("ph1.pt/mass",&pho1PtByMass);         
  diPhotonMVA->AddVariable("ph2.pt/mass",&pho2PtByMass);         
  diPhotonMVA->AddVariable("ph1.eta",&pho1Eta);             
  diPhotonMVA->AddVariable("ph2.eta",&pho2Eta);             
  diPhotonMVA->AddVariable("TMath::Cos(ph1.phi-ph2.phi)",&cosDPhi);            
  diPhotonMVA->AddVariable("ph1.idmva",&pho1IdMVA);
  diPhotonMVA->AddVariable("ph2.idmva",&pho2IdMVA);

  //book MVAs:
  diPhotonMVA->BookMVA(  methodName_diPho, weightFile_diPho);
}

void HggSelector::setBranchAddresses(){
  if(!valid) return;
  //Event info
  fChain->SetBranchAddress("lumiBlock",&lumiBlock);
  fChain->SetBranchAddress("runNumber",&runNumber);
  fChain->SetBranchAddress("evtNumber",&evtNumber);
  //fChain->SetBranchAddress("isRealData",&_isData);
<<<<<<< HEAD
  
=======
 
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
  fChain->SetBranchAddress("eeBadScFilterFlag",&eeBadScFilterFlag);
  fChain->SetBranchAddress("hcalLaserEventFilterFlag",&hcalLaserEventFilterFlag);
  fChain->SetBranchAddress("HBHENoiseFilterResultFlag",&HBHENoiseFilterResultFlag);
  fChain->SetBranchAddress("isNotDeadEcalCluster",&isNotDeadEcalCluster);
  fChain->SetBranchAddress("trackerFailureFilterFlag",&trackerFailureFilterFlag);
  fChain->SetBranchAddress("CSCHaloFilterFlag",&CSCHaloFilterFlag);
  fChain->SetBranchAddress("drDead",&drDead);
  fChain->SetBranchAddress("drBoundary",&drBoundary);
  fChain->SetBranchAddress("ECALTPFilterFlag",&ECALTPFilterFlag);

 ///information for the vertex
  fChain->SetBranchAddress("nVtx",&nVtx);
  fChain->SetBranchAddress("vtxX",vtxX);
  fChain->SetBranchAddress("vtxY",vtxY);
  fChain->SetBranchAddress("vtxZ",vtxZ);
  fChain->SetBranchAddress("vtxChi2",vtxChi2);
  fChain->SetBranchAddress("vtxNdof",vtxNdof);
  fChain->SetBranchAddress("vtxNormalizedChi2",vtxNormalizedChi2);
  fChain->SetBranchAddress("vtxTrackSize",vtxTrackSize);
  fChain->SetBranchAddress("vtxIsFake",vtxIsFake);
  fChain->SetBranchAddress("vtxIsValid",vtxIsValid);
  
  fChain->SetBranchAddress("rho", &rho);
  fChain->SetBranchAddress("rhoEtaMax44", &rhoEtaMax44);

  fChain->SetBranchAddress("pileupWeight", &pileupWeight);
 
 //objects
 fChain->SetBranchAddress("nPho",&nPho_);
 fChain->SetBranchAddress("Photons",&Photons_);

 fChain->SetBranchAddress("photonMatchedElectron",photonMatchedElectron);
 fChain->SetBranchAddress("nPair",&nPair_); 
 fChain->SetBranchAddress("ggVerticesPhotonIndices",&ggVerticesPhotonIndices);
 fChain->SetBranchAddress("ggVerticesVertexIndex",&ggVerticesVertexIndex);
 fChain->SetBranchAddress("ggVerticesPerEvtMVA",&ggVerticesPerEvtMVA);

 fChain->SetBranchAddress("nMu",&nMu_);
 fChain->SetBranchAddress("Muons",&Muons_);

 fChain->SetBranchAddress("nJet",&nJet_);
 fChain->SetBranchAddress("Jets",&Jets_);

 fChain->SetBranchAddress("nGenHiggs",&nGenHiggs);
 fChain->SetBranchAddress("GenHiggs",&GenHiggs);

 fChain->SetBranchAddress("nGenPho",&nGenPho);
 fChain->SetBranchAddress("GenPhotons",&GenPhotons);

 fChain->SetBranchAddress("nPU",&inPU);

 fChain->SetBranchAddress("PFMET",&pfMet);
 fChain->SetBranchAddress("PFMETPhi",&pfMetPhi);

<<<<<<< HEAD
=======
 fChain->SetBranchAddress("CaloMETPhi", &CaloMETPhi);
 fChain->SetBranchAddress("CaloMET", &CaloMET);

>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
 vector<string>::const_iterator trigIt;
 int i=0;
 for(trigIt=triggers.begin();trigIt!=triggers.end();trigIt++,i++){
   fChain->SetBranchAddress(trigIt->c_str(),&(triggerDec[i]));
 }

}

void HggSelector::setupOutputTree(){
  outTree = new TTree("HggOutput","");
  outTree->Branch("trigger",&trigger_,"trigger/I");
  outTree->Branch("mPair",&mPair_,"mPair/F");
  outTree->Branch("mPairNoCorr",&mPairNoCorr_,"mPairNoCorr/F");
  outTree->Branch("mPairRes",&mPairRes_,"mPairRes/F");
  outTree->Branch("mPairResWrongVtx",&mPairResWrongVtx_,"mPairResWrongVtx/F");
  outTree->Branch("diPhotonMVA",&diPhoMVA_,"diPhotonMVA/F");
  outTree->Branch("diPhotonMVAShift",diPhoMVAShift_,"diPhotonMVAShift[9]/F");
  outTree->Branch("diPhotonVtx",&diPhoVtx_,"diPhotonVtx/I");
  outTree->Branch("diPhotonVtxX",&diPhoVtxX_,"diPhotonVtxX/F");
  outTree->Branch("diPhotonVtxY",&diPhoVtxY_,"diPhotonVtxY/F");
  outTree->Branch("diPhotonVtxZ",&diPhoVtxZ_,"diPhotonVtxZ/F");
  outTree->Branch("vtxProb",&vtxProb_,"vtxProb/F");
  outTree->Branch("Mjj",&Mjj_,"Mjj");
  outTree->Branch("ptJet1",&ptJet1_,"ptJet1");
  outTree->Branch("ptJet2",&ptJet2_,"ptJet2");
  outTree->Branch("cosThetaLead",&cosThetaLead,"cosThetaLead/F");
  outTree->Branch("cat",&cat_,"cat/I");

  outTree->Branch("mPairPFCiC",&mPairPFCiC_,"mPairPFCiC/F");
  outTree->Branch("mPairNoCorrPFCiC",&mPairNoCorrPFCiC_,"mPairNoCorrPFCiC/F");
  outTree->Branch("mPairResPFCiC",&mPairResPFCiC_,"mPairResPFCiC/F");
  outTree->Branch("mPairResWrongVtxPFCiC",&mPairResWrongVtxPFCiC_,"mPairResWrongVtxPFCiC/F");
  outTree->Branch("diPhotonVtxPFCiC",&diPhoVtxPFCiC_,"diPhotonVtxPFCiC/I");
  outTree->Branch("diPhotonVtxXPFCiC",&diPhoVtxXPFCiC_,"diPhotonVtxXPFCiC/F");
  outTree->Branch("diPhotonVtxYPFCiC",&diPhoVtxYPFCiC_,"diPhotonVtxYPFCiC/F");
  outTree->Branch("diPhotonVtxZPFCiC",&diPhoVtxZPFCiC_,"diPhotonVtxZPFCiC/F");
  outTree->Branch("vtxProbPFCiC",&vtxProbPFCiC_,"vtxProbPFCiC/F");
  outTree->Branch("MjjPFCiC",&MjjPFCiC_,"MjjPFCiC");
  outTree->Branch("ptJet1PFCiC",&ptJet1PFCiC_,"ptJet1PFCiC");
  outTree->Branch("ptJet2PFCiC",&ptJet2PFCiC_,"ptJet2PFCiC");
  outTree->Branch("cosThetaLeadPFCiC",&cosThetaLeadPFCiC,"cosThetaLeadPFCiC/F");
  outTree->Branch("catPFCiC",&catPFCiC_,"catPFCiC/I");

  outTree->Branch("eeBadScFilterFlag",&eeBadScFilterFlagOut);
  outTree->Branch("hcalLaserEventFilterFlag",&hcalLaserEventFilterFlagOut);
  outTree->Branch("HBHENoiseFilterResultFlag",&HBHENoiseFilterResultFlagOut);
  outTree->Branch("isNotDeadEcalCluster",&isNotDeadEcalClusterOut);
  outTree->Branch("trackerFailureFilterFlag",&trackerFailureFilterFlagOut);
  outTree->Branch("CSCHaloFilterFlag",&CSCHaloFilterFlagOut);
  outTree->Branch("drDead",&drDeadOut);
  outTree->Branch("drBoundary",&drBoundaryOut);
  outTree->Branch("ECALTPFilterFlag",&ECALTPFilterFlagOut);

  outTree->Branch("mPairCiC",&mPairCiC_,"mPairCiC/F");
  outTree->Branch("mPairNoCorrCiC",&mPairNoCorrCiC_,"mPairNoCorrCiC/F");
  outTree->Branch("mPairResCiC",&mPairResCiC_,"mPairResCiC/F");
  outTree->Branch("mPairResWrongVtxCiC",&mPairResWrongVtxCiC_,"mPairResWrongVtxCiC/F");
  outTree->Branch("diPhotonVtxCiC",&diPhoVtxCiC_,"diPhotonVtxCiC/I");
  outTree->Branch("diPhotonVtxXCiC",&diPhoVtxXCiC_,"diPhotonVtxXCiC/F");
  outTree->Branch("diPhotonVtxYCiC",&diPhoVtxYCiC_,"diPhotonVtxYCiC/F");
  outTree->Branch("diPhotonVtxZCiC",&diPhoVtxZCiC_,"diPhotonVtxZCiC/F");
  outTree->Branch("vtxProbCiC",&vtxProbCiC_,"vtxProbCiC/F");
  outTree->Branch("MjjCiC",&MjjCiC_,"MjjCiC");
  outTree->Branch("ptJet1CiC",&ptJet1CiC_,"ptJet1CiC");
  outTree->Branch("ptJet2CiC",&ptJet2CiC_,"ptJet2CiC");
  outTree->Branch("cosThetaLeadCiC",&cosThetaLeadCiC,"cosThetaLeadCiC/F");

  outTree->Branch("nPhoton",&nOutPhotons_,"nPhoton/I");
  outTree->Branch("Photon",&OutPhotons_);
  outTree->Branch("nPhotonPFCiC",&nOutPhotonsPFCiC_,"nPhotonPFCiC/I");
  outTree->Branch("PhotonPFCiC",&OutPhotonsPFCiC_);
  outTree->Branch("nPhotonCiC",&nOutPhotonsCiC_,"nPhotonCiC/I");
  outTree->Branch("PhotonCiC",&OutPhotonsCiC_);

  outTree->Branch("MET",&MET);
  outTree->Branch("METPhi",&METPhi);

  outTree->Branch("genHiggsPt",&genHiggsPt);
  outTree->Branch("genHiggsVx",&genHiggsVx);
  outTree->Branch("genHiggsVy",&genHiggsVy);
  outTree->Branch("genHiggsVz",&genHiggsVz);
  
  outTree->Branch("evtWeight",&evtWeight);

  outTree->Branch("ptGenPho1",&ptGenPho1);
  outTree->Branch("etaGenPho1",&etaGenPho1);
  outTree->Branch("phiGenPho1",&phiGenPho1);
  outTree->Branch("energyGenPho1",&energyGenPho1);

  outTree->Branch("ptGenPho2",&ptGenPho2);
  outTree->Branch("etaGenPho2",&etaGenPho2);
  outTree->Branch("phiGenPho2",&phiGenPho2);
  outTree->Branch("energyGenPho2",&energyGenPho2);

  outTree->Branch("mPairScale",&mPairScale);
  outTree->Branch("pho1MVAScale",&pho1MVAScale);
  outTree->Branch("pho2MVAScale",&pho2MVAScale);
  outTree->Branch("diPhoMVAScale",&diPhoMVAScale);
  outTree->Branch("mPairSmear",&mPairSmear);
  outTree->Branch("pho1MVASmear",&pho1MVASmear);
  outTree->Branch("pho2MVASmear",&pho2MVASmear);
  outTree->Branch("diPhoMVASmear",&diPhoMVASmear);

  outTree->Branch("mPairScalePFCiC",&mPairScalePFCiC);
  outTree->Branch("mPairSmearPFCiC",&mPairSmearPFCiC);

  outTree->Branch("mPairScaleCiC",&mPairScaleCiC);
  outTree->Branch("mPairSmearCiC",&mPairSmearCiC);

  outTree->Branch("nPU",&nPU_);
  outTree->Branch("nVtx",&nVtxOut);

  outTree->Branch("runNumber",&runNumberOut);
  outTree->Branch("evtNumber",&evtNumberOut);
  outTree->Branch("lumiBlock",&lumiBlockOut);
  

  if(doMuMuGamma){
    outTreeMuMuG = new TTree("MuMuGamma","");
    outTreeMuMuG->Branch("nMuMuG",&nMuMuG);

    outTreeMuMuG->Branch("massMuMuGamma",massMuMuGamma,"massMuMuGamma[nMuMuG]");
    outTreeMuMuG->Branch("massMuMuRegGamma",massMuMuRegGamma,"massMuMuRegGamma[nMuMuG]");
    outTreeMuMuG->Branch("massMuMuScaleGamma",massMuMuScaleGamma,"massMuMuScaleGamma[nMuMuG]");
    outTreeMuMuG->Branch("massMuMuGenGamma",massMuMuGenGamma,"massMuMuGenGamma[nMuMuG]");
    outTreeMuMuG->Branch("massMuMu",massMuMu,"massMuMu[nMuMuG]");
    outTreeMuMuG->Branch("puWeight",puWeight,"puWeight[nMuMuG]");

    outTreeMuMuG->Branch("Muon1",&MMG_Mu1);
    outTreeMuMuG->Branch("Muon2",&MMG_Mu2);
    outTreeMuMuG->Branch("Photon",&MMG_Pho);
  
    outTreeMuMuG->Branch("isosumoetPho",isosumoetPho,"isosumoetPho[nMuMuG]");
    outTreeMuMuG->Branch("mvaPho",mvaPho,"mvaPho[nMuMuG]");
  }
<<<<<<< HEAD
=======

  ////////////RAZOR VARIABLES////////////////
  //razor
  outTree->Branch("nBtags",&nBtags,"nBtags/I");
  outTree->Branch("nJets",&nJets,"nJets/I");
  outTree->Branch("bothBarrel", &bothBarrel, "bothBarrel/I");

  outTree->Branch("CaloMETPhi", &CaloMETPhi, "CaloMETPhi/F");
  outTree->Branch("CaloMET", &CaloMET, "CaloMET/F");

  //NO SEEDING VARIABLES
  outTree->Branch("PFMR",&PFMR,"PFMR/F");
  outTree->Branch("PFR",&PFR,"PFR/F");

  outTree->Branch("ptHem1",&ptHem1,"ptHem1/F");
  outTree->Branch("etaHem1",&etaHem1,"etaHem1/F");
  outTree->Branch("phiHem1",&phiHem1,"phiHem1/F");
  
  outTree->Branch("ptHem2",&ptHem2,"ptHem2/F");
  outTree->Branch("etaHem2",&etaHem2,"etaHem2/F");
  outTree->Branch("phiHem2",&phiHem2,"phiHem2/F");  

  outTree->Branch("mHem1",&mHem1,"mHem1/F");
  outTree->Branch("mHem2",&mHem2,"mHem2/F");

  //SAME SIDE VARIABLES
  outTree->Branch("PFMR_SS",&PFMR_SS,"PFMR_SS/F");
  outTree->Branch("PFR_SS",&PFR_SS,"PFR_SS/F");

  outTree->Branch("ptHem1_SS",&ptHem1_SS,"ptHem1_SS/F");
  outTree->Branch("etaHem1_SS",&etaHem1_SS,"etaHem1_SS/F");
  outTree->Branch("phiHem1_SS",&phiHem1_SS,"phiHem1_SS/F");
  
  outTree->Branch("ptHem2_SS",&ptHem2_SS,"ptHem2_SS/F");
  outTree->Branch("etaHem2_SS",&etaHem2_SS,"etaHem2_SS/F");
  outTree->Branch("phiHem2_SS",&phiHem2_SS,"phiHem2_SS/F");  
  
  outTree->Branch("mHem1_SS",&mHem1_SS,"mHem1_SS/F");
  outTree->Branch("mHem2_SS",&mHem2_SS,"mHem2_SS/F");

  //OPPOSITE SIDE VARIABLES
  outTree->Branch("PFMR_OS",&PFMR_OS,"PFMR_OS/F");
  outTree->Branch("PFR_OS",&PFR_OS,"PFR_OS/F");

  outTree->Branch("ptHem1_OS",&ptHem1_OS,"ptHem1_OS/F");
  outTree->Branch("etaHem1_OS",&etaHem1_OS,"etaHem1_OS/F");
  outTree->Branch("phiHem1_OS",&phiHem1_OS,"phiHem1_OS/F");
  
  outTree->Branch("ptHem2_OS",&ptHem2_OS,"ptHem2_OS/F");
  outTree->Branch("etaHem2_OS",&etaHem2_OS,"etaHem2_OS/F");
  outTree->Branch("phiHem2_OS",&phiHem2_OS,"phiHem2_OS/F");  

  outTree->Branch("mHem1_OS",&mHem1_OS,"mHem1_OS/F");
  outTree->Branch("mHem2_OS",&mHem2_OS,"mHem2_OS/F");
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
}

void HggSelector::fillMuMuGamma(){
  nMuMuG=0;
  TVector3 vtx(vtxX[0],vtxY[0],vtxZ[0]);
  for(int iMu1=0;iMu1<nMu_;iMu1++){
    VecbosMu mu1 = Muons_->at(iMu1);    
    if(mu1.pt < 10) continue;
    if(!mu1.isGlobalMuon || !mu1.isTrackerMuon) continue;
    if(mu1.nTrackHits <= 10 || mu1.nPixelHits==0) continue;
    if(mu1.trackImpactPar >=0.2) continue;
    if(mu1.trkIso >= 3) continue;
    for(int iMu2=iMu1+1;iMu2<nMu_;iMu2++){
      VecbosMu mu2 = Muons_->at(iMu2);
      if(mu2.pt < 10) continue;
      if(mu1.pt<20 && mu2.pt<20) continue; // 20/10 selection
      if(mu1.charge*mu2.charge >=0) continue; //opposite charge
      if(!mu2.isGlobalMuon || !mu2.isTrackerMuon) continue;
      if(mu2.nTrackHits <= 10 || mu2.nPixelHits==0) continue;
      if(mu2.trackImpactPar >=0.2) continue;
      if(mu2.trkIso >= 3) continue;
      for(int iPho=0; iPho<nPho_;iPho++){
	VecbosPho pho = Photons_->at(iPho);

	//fill the different masses
	TLorentzVector p4Mu1; p4Mu1.SetPtEtaPhiM(mu1.pt,mu1.eta,mu1.phi,0.106);
	TLorentzVector p4Mu2; p4Mu2.SetPtEtaPhiM(mu2.pt,mu2.eta,mu2.phi,0.106);

	massMuMu[nMuMuG] = (p4Mu1+p4Mu2).M();
	massMuMuGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.energy,false)).M();
	massMuMuRegGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.correctedEnergy,false)).M();
	massMuMuScaleGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.scaledEnergy,false)).M();
	if(pho.genMatch.index>=0) massMuMuGenGamma[nMuMuG] = (p4Mu1+p4Mu2+pho.p4FromVtx(vtx,pho.genMatch.energy,false)).M();
	else massMuMuGenGamma[nMuMuG] = -1;

	MMG_Mu1.push_back(mu1);
	MMG_Mu2.push_back(mu2);
	MMG_Pho.push_back(pho);
	float eT = pho.p4FromVtx(vtx,pho.energy,false).Et();
	isosumoetPho[nMuMuG] = (pho.dr03EcalRecHitSumEtCone + pho.dr04HcalTowerSumEtCone + pho.dr03TrkSumPtHollowCone + isoSumConst - rho*rhoFac)/eT;
	mvaPho[nMuMuG] = PhotonID->getIdMVA(&pho,nVtx,rho,0); 
	
	nMuMuG++;
      }
    }
  }


  if(doMuMuGamma) outTreeMuMuG->Fill();
}

#define debugMjj 0
float HggSelector::getVBFMjj(VecbosPho* pho1, VecbosPho* pho2,TVector3 SelVtx,float *jetPts){
  TLorentzVector p1 = pho1->p4FromVtx(SelVtx,pho1->finalEnergy);
  TLorentzVector p2 = pho2->p4FromVtx(SelVtx,pho2->finalEnergy);
  jetPts[0]=0.; jetPts[1]=0.;
  if( max(p1.Pt(),p2.Pt())/(p1+p2).M() < 60./120. ) return 0;
  if( min(p1.Pt(),p2.Pt()) < 25. ) return 0;
  //if( (p1+p2).M() < 120. ) return 0; 
  if(debugMjj) cout << "NJets: " << nJet_ << endl;
  if(nJet_<2) return 0;
  int i1=-1,i2=-1;
  float maxSumPt=0;
  VecbosJet *jet1=0, *jet2=0;
  std::vector<VecbosJet>::iterator j1It;
  std::vector<VecbosJet>::iterator j2It;
  for(j1It = Jets_->begin(); j1It != Jets_->end(); j1It++){
    if(j1It->pt < 20.) continue;
    if(!passJetID(&*j1It)) continue;
    if(DeltaR(j1It->eta,p1.Eta(),j1It->phi,p1.Phi()) < 0.5 ) continue;
    if(DeltaR(j1It->eta,p2.Eta(),j1It->phi,p2.Phi()) < 0.5 ) continue;
    for(j2It = j1It+1; j2It != Jets_->end(); j2It++){
      if(j2It->pt < 20.) continue;
      if(!passJetID(&*j2It)) continue;
      if(DeltaR(j2It->eta,p1.Eta(),j2It->phi,p1.Phi()) < 0.5 ) continue;
      if(DeltaR(j2It->eta,p2.Eta(),j2It->phi,p2.Phi()) < 0.5 ) continue;

      if(j1It->pt < 30. && j2It->pt < 30.) continue;
      if(j1It->pt + j2It->pt > maxSumPt){
	maxSumPt = j1It->pt + j2It->pt;
	jet1 = &*j1It; jet2 = &*j2It;
      }
    }
  }

  if(jet1==0 || jet2==0) return 0;
  TLorentzVector j1 = jet1->getP4();
  TLorentzVector j2 = jet2->getP4();
  MjjDists["BeforeVBF"]->Fill((j1+j2).M());
  float dEtaJ = fabs(j1.Eta()-j2.Eta());
  if(debugMjj) cout << "dEtaJ: " << dEtaJ << endl;
  MjjDists["dEta"]->Fill(dEtaJ);
  if(dEtaJ < 3.) return 0;
  MjjDists["AfterDEta"]->Fill( (j1+j2).M() );
  TLorentzVector ggSystem = p1+p2;
  float Z = ggSystem.Eta() - (j1.Eta()+j2.Eta())/2.;
  if(debugMjj) cout << "Z: " << Z << endl;
  MjjDists["Z"]->Fill(Z);
  if(fabs(Z)>2.5) return 0;
  MjjDists["AfterZ"]->Fill( (j1+j2).M() );
  TLorentzVector jjSystem = j1+j2;
  if(debugMjj) cout << "dPhi jj gg: " << fabs(jjSystem.DeltaPhi(ggSystem)) << endl;
  MjjDists["dPhi_jj_gg"]->Fill( fabs(jjSystem.DeltaPhi(ggSystem)) );
  if( fabs(jjSystem.DeltaPhi(ggSystem)) < 2.6 ) return 0;
  MjjDists["Final"]->Fill( (j1+j2).M() );
  jetPts[0] = j1.Pt();
  jetPts[1] = j2.Pt();
  if(debugMjj) cout << "jj Mass: " << jjSystem.M() <<endl;
  return jjSystem.M();
}

bool HggSelector::passJetID(VecbosJet* jet){
  const int nJetCat=4;
  const float maxJetEta[nJetCat] = {2.5,2.75,3,4.7};
  const float betaStarSlope[nJetCat] = {0.2,0.3,999,999};
  const float rmsCut[nJetCat] = {0.06,0.05,0.05,0.055};

  int JetCat=-1;
  for(int i=0;i<nJetCat;i++){
    if(jet->eta < maxJetEta[i]){
      JetCat = i; break;
    }
  }
  if(JetCat == -1) return false;

  if(jet->betaStarClassicIdMVA > betaStarSlope[JetCat]*TMath::Log(nVtx-0.64)) return false; //warning, I think this is WRONG! but its how MIT does it...
  if(jet->rmsCandsHand > rmsCut[JetCat]) return false;
  return true;
}

int HggSelector::getCategory(){
  if(diPhoMVA_ < -0.05) return -1;
  if(Mjj_ >=500 && ptJet1_ >= 30. && ptJet2_ >= 30.) return 4;
  if(Mjj_ >=250 && (ptJet1_ >= 30. || ptJet2_ >= 30.)) return 5;
  
  if(diPhoMVA_ > 0.88) return 0;
  if(diPhoMVA_ > 0.71) return 1;
  if(diPhoMVA_ > 0.50) return 2;
  return 3;
}

int HggSelector::getCategoryPFCiC(){
  if(Mjj_ >=500 && ptJet1_ >= 30. && ptJet2_ >= 30.) return 4;
  if(Mjj_ >=250 && (ptJet1_ >= 30. || ptJet2_ >= 30.)) return 5;

  bool EBEB = (fabs(pho1_.SC.eta) < 1.48) && (fabs(pho2_.SC.eta) < 1.48);
  bool R9R9 = (pho1_.SC.r9 > 0.94) && (pho2_.SC.r9 > 0.94);
  
  if(R9R9 && EBEB) return 0;
  if(!R9R9 && EBEB) return 1;
  if(R9R9 && !EBEB) return 2;
  if(!R9R9 && !EBEB) return 3;
  return -1;
}

bool HggSelector::requireTrigger(){
  if(!_isData) return true; //no triggers on MC

  if(triggers.size()==0) return true; //no trigger selection
  
  for(int i=0;i<triggers.size();i++){
    if(triggerDec[i]) return true;
  }
  return false;
}


int HggSelector::getVertexIndex(int indexPho1,int indexPho2){
  if(forceVtxZero) return 0;
  //get the vertex selected by the di photon vertex MVA
  int selectedVertex = -1;
  if(indexPho1 > nPho_ || indexPho2 > nPho_) return -1;
  int origIndex1 = Photons_->at(indexPho1).index;  //index is based on the original photon index
  int origIndex2 = Photons_->at(indexPho2).index;  
  if(debugSelector) cout << "Original Indices:" << origIndex1 << "  " << origIndex2 << endl;
  for(int i=0;i<nPair_;i++){
    //the first of these should be the correct order, but just in case ....
    if(ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex1,origIndex2) ||
       ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex2,origIndex1)){
      selectedVertex = ggVerticesVertexIndex->at(i).first;
      break;
    }
  }
  if(debugSelector) cout << "Returning Vertex Index: " << selectedVertex <<endl;
  return selectedVertex;
}
float HggSelector::getVertexMVA(int indexPho1,int indexPho2){
  //get the vertex selected by the di photon vertex MVA
  int origIndex1 = Photons_->at(indexPho1).index;  //index is based on the original photon index
  int origIndex2 = Photons_->at(indexPho2).index;  
  float selectedVertex = -1;
  for(int i=0;i<nPair_;i++){
    //the first of these should be the correct order, but just in case ....
    if(ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex1,origIndex2) ||
       ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex2,origIndex1)){
      selectedVertex = ggVerticesVertexIndex->at(i).second;
      break;
    }
  }
  return selectedVertex;
}

float HggSelector::getVertexProb(int indexPho1,int indexPho2){
  int origIndex1 = Photons_->at(indexPho1).index;  //index is based on the original photon index
  int origIndex2 = Photons_->at(indexPho2).index;  
  float evtmva = -1e6;
  for(int i=0;i<nPair_;i++){
    //the first of these should be the correct order, but just in case ....
    if(ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex1,origIndex2) ||
       ggVerticesPhotonIndices->at(i) == std::pair<int,int>(origIndex2,origIndex1)){
      evtmva = ggVerticesPerEvtMVA->at(i);
      break;
    }
  }
  return vtxprob=1.-0.49*(evtmva+1.0);
  
}

float HggSelector::getDiPhoMVA(int indexPho1, int indexPho2, float mva1, float mva2, bool usePFSC){
  if(debugSelector) cout << "getDiPhoMVA" <<endl;  
  int selectedVertex = getVertexIndex(indexPho1,indexPho2);
  if(selectedVertex<0 || selectedVertex >= nVtx){
    cout << "WARNING: Photons " << indexPho1 << " and " << indexPho2 
	 << " have no selected vertex!" << endl
	 << "Skipping this pair" << endl;
    return -9999.;
  }
  if(debugSelector) cout << mva1 << "  " << mva2 << endl;
  if(mva1 < -999 || mva2 < -999) return -9999.;
  if(debugSelector) cout << "getting photons" <<endl;  
  VecbosPho pho1 = Photons_->at(indexPho1);
  VecbosPho pho2 = Photons_->at(indexPho2);
  if(debugSelector) cout << "selected Vertex: " << selectedVertex <<endl;  

  if(debugSelector){
    std::cout << "pho1: " << pho1.finalEnergy << "  " << pho1.scaledEnergy << "  " << pho1.dEoE << "  " << pho1.eta << "  " << std::endl;
    std::cout << "pho2: " << pho2.finalEnergy << "  " << pho2.scaledEnergy << "  " << pho2.dEoE << "  " << pho2.eta << "  " << std::endl;
  }
  TVector3 vtxPos(vtxX[selectedVertex],vtxY[selectedVertex],vtxZ[selectedVertex]);

  TLorentzVector p1 = pho1.p4FromVtx(vtxPos,pho1.finalEnergy);
  TLorentzVector p2 = pho2.p4FromVtx(vtxPos,pho2.finalEnergy);

  VecbosPho* phoLead = (p1.Et() > p2.Et() ? &pho1 : & pho2); //select the leading photon
  VecbosPho* phoSubLead = (p1.Et() > p2.Et() ? &pho2 : & pho1); //select the leading photon

  VecbosSC *scLead;
  VecbosSC *scSubLead;
  if(usePFSC){
    scLead = &(phoLead->PFSC);
    scSubLead = &(phoSubLead->PFSC);
  }else{
    scLead = &(phoLead->SC);
    scSubLead = &(phoSubLead->SC);
  }

  float mvaLead = (p1.Et() > p2.Et() ? mva1 : mva2);
  float mvaSubLead = (p1.Et() > p2.Et() ? mva2 : mva1);

  float mPair = (p1+p2).M();
  if(debugSelector) std::cout << "mPair: " << mPair << std::endl;
  //fill variables
  smearedMassErrByMass = massRes->getMassResolutionEonly(&pho1,&pho2,vtxPos)/mPair;
  smearedMassErrByMassWrongVtx = massRes->getMassResolution(&pho1,&pho2,vtxPos,true)/mPair;
  if(debugSelector) cout << "Mass Erro resolved. Selected Vertex: " << selectedVertex << endl;

  vtxprob= getVertexProb(indexPho1,indexPho2);


  if(debugSelector) cout << "vtx prob: " << vtxprob << endl;
  pho1PtByMass = max(p1.Et(),p2.Et())/mPair;
  pho2PtByMass = min(p1.Et(),p2.Et())/mPair;
  if(debugSelector) cout << "pho1/m: " << pho1PtByMass << endl;
  if(debugSelector) cout << "pho2/m: " << pho2PtByMass << endl;
  pho1Eta = (p1.Et() > p2.Et() ? p1.Eta(): p2.Eta());
  pho2Eta = (p1.Et() > p2.Et() ? p2.Eta(): p1.Eta());
  if(debugSelector) cout << "pho1 eta: " << pho1Eta << endl;
  if(debugSelector) cout << "pho2 eta: " << pho2Eta << endl;
  cosDPhi = TMath::Cos(DeltaPhi(scLead->phi,scSubLead->phi));
  if(debugSelector) cout << "cos dPhi: " << cosDPhi << endl;
  pho1IdMVA = mvaLead;
  pho2IdMVA = mvaSubLead;

  MVAInputs["massRes"]->Fill(smearedMassErrByMass);
  MVAInputs["massResWrongVtx"]->Fill(smearedMassErrByMassWrongVtx);
  MVAInputs["vtxProb"]->Fill(vtxprob);
  MVAInputs["p1_EtByM"]->Fill(pho1PtByMass);
  MVAInputs["p2_EtByM"]->Fill(pho2PtByMass);
  MVAInputs["p1_Eta"]->Fill(pho1Eta);
  MVAInputs["p2_Eta"]->Fill(pho2Eta);
  MVAInputs["CosDPhi"]->Fill(cosDPhi);
  MVAInputs["p1_idMVA"]->Fill(pho1IdMVA);
  MVAInputs["p2_idMVA"]->Fill(pho2IdMVA);

  return diPhotonMVA->EvaluateMVA(methodName_diPho);
}
<<<<<<< HEAD
=======
  ////////////RAZOR METHODS////////////////
vector<TLorentzVector> HggSelector::CombineJets_R_no_seed(vector<TLorentzVector> myjets,TLorentzVector ph1, TLorentzVector ph2){
  vector<TLorentzVector> mynewjets;
  TLorentzVector j1, j2;
  
  myjets.push_back(ph1);
  myjets.push_back(ph2);

  int N_comb = 1;
  for(int i = 0; i < myjets.size(); i++){
    N_comb *= 2;
  }

  double M_min = 9999999999.0;
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TLorentzVector j_temp1, j_temp2;    
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
        j_temp1 += myjets[count];
      } else {
        j_temp2 += myjets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }    
    double M_temp = j_temp1.M2()+j_temp2.M2();
    // smallest mass
    if(M_temp < M_min){
      // R selection
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }

  // set masses to 0
  j1.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),0.0);
  j2.SetPtEtaPhiM(j2.Pt(),j2.Eta(),j2.Phi(),0.0);
  
  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }
  
  mynewjets.push_back(j1);
  mynewjets.push_back(j2);
  return mynewjets;    
}

vector<TLorentzVector> HggSelector::CombineJets_R_SSorOS(vector<TLorentzVector> myjets,TLorentzVector ph1, TLorentzVector ph2, bool SS){
  
  vector<TLorentzVector> mynewjets;
  TLorentzVector j1, j2;
  
  int N_comb = 1;
  for(int i = 0; i < myjets.size(); i++){
    N_comb *= 2;
  }

  double M_min = 9999999999.0;
  int j_count;
  for(int i = 0; i < N_comb; i++){
    TLorentzVector j_temp1, j_temp2;

    //seed the photons based on SS
    if(SS) {
      j_temp1+=ph1;
      j_temp1+=ph2;
    }
    else {
      j_temp1+=ph1;
      j_temp2+=ph2;
    }

    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
        j_temp1 += myjets[count];
      } else {
        j_temp2 += myjets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }    
    double M_temp = j_temp1.M2()+j_temp2.M2();
    // smallest mass
    if(M_temp < M_min){
      // R selection
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }

  //handle the special case of jetlist 1
  if(myjets.size() == 1 && !SS) {
    float poss1 = (ph1 + myjets[0]).M2() + ph2.M2();
    float poss2 = ph1.M2() + (ph2 + myjets[0]).M2();
    
    if(poss1 < poss2) {
      j1 = ph1 + myjets[0];
      j2 = ph2;
    }
    else {
      j1 = ph2 + myjets[0];
      j2 = ph1;
    }
  }
  if(myjets.size() == 1 && SS) {
    j1 = ph1 + ph2;
    j2 = myjets[0];
  }
  
  //if we require no jet we end up with terrible SS scenarios
  if(myjets.size() == 0) {
    if(SS) {
      cout << "WARNING: NO JETS AND SAME SAME MR Rsq RECONSRUCTION" << endl;
      j1 = ph1 + ph2;
      j2.SetPtEtaPhiM(0, 0, 0, 0.0);;
    }
    else {
      j1 = ph1;
      j2 = ph2;
    }
  }

  // set masses to 0
  j1.SetPtEtaPhiM(j1.Pt(),j1.Eta(),j1.Phi(),0.0);
  j2.SetPtEtaPhiM(j2.Pt(),j2.Eta(),j2.Phi(),0.0);
  

  if(j2.Pt() > j1.Pt()){
    TLorentzVector temp = j1;
    j1 = j2;
    j2 = temp;
  }
  
  mynewjets.push_back(j1);
  mynewjets.push_back(j2);
  return mynewjets;  
}

double HggSelector::CalcGammaMRstar(TLorentzVector ja, TLorentzVector jb){
  double A = ja.P();
  double B = jb.P();
  double az = ja.Pz();
  double bz = jb.Pz();
  TVector3 jaT, jbT;
  jaT.SetXYZ(ja.Px(),ja.Py(),0.0);
  jbT.SetXYZ(jb.Px(),jb.Py(),0.0);
  double ATBT = (jaT+jbT).Mag2();

  double temp = sqrt((A+B)*(A+B)-(az+bz)*(az+bz)-
                     (jbT.Dot(jbT)-jaT.Dot(jaT))*(jbT.Dot(jbT)-jaT.Dot(jaT))/(jaT+jbT).Mag2());

  double mybeta = (jbT.Dot(jbT)-jaT.Dot(jaT))/
    sqrt(ATBT*((A+B)*(A+B)-(az+bz)*(az+bz)));

  double mygamma = 1./sqrt(1.-mybeta*mybeta);

  //gamma times MRstar                                                                                                                                                                              
  temp *= mygamma;

  return temp;
}

double HggSelector::CalcMTR(TLorentzVector ja, TLorentzVector jb, TVector3 met){

  double temp = met.Mag()*(ja.Pt()+jb.Pt()) - met.Dot(ja.Vect()+jb.Vect());
  temp /= 2.;

  temp = sqrt(temp);

  return temp;
}

vector<TLorentzVector> HggSelector::GetJetList(TLorentzVector p1, TLorentzVector p2){
  vector<TLorentzVector> jetlist;

  std::vector<VecbosJet>::iterator jIt;
  for(jIt = Jets_->begin(); jIt != Jets_->end(); jIt++){
    //eta cut
    if(fabs(jIt->eta) > 2.6) continue;
    //pt cut
    if(jIt->pt < 30.) continue;
    //jet id cut
    if(!passJetID(&*jIt)) continue;

    //photon matching
    if(DeltaR(jIt->eta,p1.Eta(),jIt->phi,p1.Phi()) < 0.5 ) continue;
    if(DeltaR(jIt->eta,p2.Eta(),jIt->phi,p2.Phi()) < 0.5 ) continue;

    jetlist.push_back(jIt->getP4());
  }

  return jetlist;
}

void HggSelector::PrintEventNumbers() {
  cout << runNumber << ":" << lumiBlock << ":" << evtNumber << endl;
}

bool HggSelector::PassMETFilters(){
  //only using MET filters Javier is using (bits 0 3 4 6 7 8 respectively

  // remove bits 6 7 8 for now (unfilled) 
  bool decision =   (ECALTPFilterFlag && CSCHaloFilterFlag && trackerFailureFilterFlag && HBHENoiseFilterResultFlag && hcalLaserEventFilterFlag && eeBadScFilterFlag);  
  //bool decision =   (ECALTPFilterFlag && CSCHaloFilterFlag && trackerFailureFilterFlag);
  if( !decision ) {
    cout << "------------Begin MET FLAG-----------" << endl;
  }

  //print flags which occur
  if( !ECALTPFilterFlag ) cout << "Ecal dead cell Flagged Bit 0";PrintEventNumbers();
  if( !CSCHaloFilterFlag ) cout << "CSC Beam Halo Flagged Bit 3";PrintEventNumbers();
  if( !trackerFailureFilterFlag ) cout << "tracker Failure Flagged Bit 4";PrintEventNumbers();
  if( !HBHENoiseFilterResultFlag ) cout << "HBHE Noise Flagged Bit 6" ;PrintEventNumbers();
  if( !hcalLaserEventFilterFlag ) cout << "HCAL Laser Flagged Bit 7" ;PrintEventNumbers();
  if( !eeBadScFilterFlag ) cout << "EE Bad SC Flagged Bit 8" ;PrintEventNumbers();
  if( !decision ) cout << "--------------End MET FLAG-------------" << endl;
  return decision;
}

void HggSelector::FillRazorVarsWith(int n){
  nBtags = 0;        
  
  //no seed
  PFMR = n;
  PFR = n;
  ptHem1= n ;
  etaHem1= n ;
  phiHem1= n ;        
  ptHem2= n ;
  etaHem2= n ;
  phiHem2= n ;      
  mHem1=n;
  mHem2=n;

  //same side
  PFMR_SS = n;
  PFR_SS = n;
  ptHem1_SS= n ;
  etaHem1_SS= n ;
  phiHem1_SS= n ;        
  ptHem2_SS= n ;
  etaHem2_SS= n ;
  phiHem2_SS= n ;      
  mHem1_SS=n;
  mHem2_SS=n;
  //opposite side
  PFMR_OS = n;
  PFR_OS = n;
  ptHem1_OS= n ;
  etaHem1_OS= n ;
  phiHem1_OS= n ;        
  ptHem2_OS= n ;
  etaHem2_OS= n ;
  phiHem2_OS= n ;      
  mHem1_OS=n;
  mHem2_OS=n;
}

//after the initeventflag is called we can check the bad event
//list with this method
bool HggSelector::isFlagged(){
  EventIndex index;
  
  //fixed the names of variables
  index.EventNumber = evtNumber;
  index.RunNumber = runNumber;
  
  if(EventCounts.find(index) == EventCounts.end()){ //yes
    //cout << "Event not in list" << endl;
    return false;
  } else {
    //cout << "Event IS in list" << endl;
    return true;
    if(EventCounts[index] == 1){
      //cout << "Found the event - all is good in the world" << endl;
      return true;
    }
  } 
}

//initialize the list of bad events
//and construct the eventindex object
void HggSelector::InitEventFlag(char *s_Event){
  ifstream inputfile(s_Event);
  
  int RUN_NUMBER;
  int LS_NUMBER;
  Long64_t EVENT_NUMBER;
  
  EventIndex index;
  
  cout << "Reading bad event list" << endl;
	
  while(!inputfile.eof()){
    inputfile >> RUN_NUMBER >> LS_NUMBER >> EVENT_NUMBER;
    
    index.RunNumber = RUN_NUMBER;
    index.EventNumber = EVENT_NUMBER;
    
    if(index.RunNumber < 0 || index.EventNumber < 0)
      continue;
    
    //Is this event/run-number combo alreading in the map?
    if(EventCounts.find(index) == EventCounts.end()){ //no
      EventCounts.insert(pair<EventIndex, int>(index, 1));
    } 
    else { //yes
      EventCounts[index] = EventCounts[index] + 1;
    }
    
    index.RunNumber = -1;
    index.EventNumber = -1;
  }
}

//The event index object for checking the list of bad

  
>>>>>>> d9db9bdd2e20dd5561b4790426fe877a8407754a
