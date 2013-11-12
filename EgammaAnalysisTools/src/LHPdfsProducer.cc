#include <iostream>
#include <string> 
#include <math.h>

#include "TLorentzVector.h"

#include "CommonTools/include/Utils.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/RedEleIDTree.hh"
#include "EgammaAnalysisTools/include/LHPdfsProducer.hh"

#include "cajun/json/reader.h"
#include "cajun/json/elements.h"
#include <CommonTools/include/TriggerMask.hh>

using namespace bits;
using namespace std;

LHPdfsProducer::LHPdfsProducer(TTree *tree)
  : EgammaBase(tree) {

  jsonFile = "";
  lastFile = "";
  
  std::string fileCuts("config/LHPdfsProducer/cuts.txt");
  std::string fileSwitches("config/LHPdfsProducer/switches.txt");
  
  m_selection = new Selection(fileCuts,fileSwitches);
  m_counters  = new Counters();
  configSelection(m_selection,m_counters);

  // data or MC
  isData_ = m_selection->getSwitch("isData");

  // reading GoodRUN LS 
  std::cout << "[GoodRunLS]::goodRunLS is " << m_selection->getSwitch("goodRunLS") << " isData is " <<  m_selection->getSwitch("isData") << std::endl;

  // to read good run list
  if (m_selection->getSwitch("goodRunLS") && m_selection->getSwitch("isData")) {
    std::string goodRunGiasoneFile = "config/LHPdfsProducer/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }

  // single electron efficiency
  TString selectionString(m_selection->getStringParameter("electronIDType"));
  cout << "=== CONFIGURING " << selectionString << " SYMMETRIC ELECTRON ID ===" << endl;
  EgammaCutBasedID.ConfigureNoClass("config/"+selectionString);
  EgammaCutBasedID.ConfigureEcalCleaner("config/");
}

LHPdfsProducer::~LHPdfsProducer() { 
  
  delete m_selection;
  delete m_counters;
}


// PDF for probe electrons within the acceptance and loose isolated in the tracker
// the tag is the one with the best match to the Z
// the tag must be within the acceptance, tracker isolated and loose identified
void LHPdfsProducer::LoopZTagAndProbe(const char *treefilesuffix) {
  
  if(fChain == 0) return;
  
  bookHistos();

  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesSignal();
  reducedTree.addCategories();
  reducedTree.addMore();        // to find the best cut
  
  // counters
  int allevents   = 0;
  int trigger     = 0;
  int twoele      = 0;
  int eleTot      = 0;
  int eleEta      = 0;
  int elePt       = 0;
  int invmass     = 0;
  int tagId       = 0;
  int tagIsol     = 0;
  int probeIsol   = 0;
  int tagProbeTot = 0;

  // json 
  unsigned int lastLumi=0;
  unsigned int lastRun=0;

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();

  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // reload trigger mask
    bool newTriggerMask = false;
    if(isData_) newTriggerMask = true;
    reloadTriggerMask(newTriggerMask);
    
    // Good Run selection 
    if (isData_ && m_selection->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun  = runNumber;
        lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (isData_ && m_selection->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    allevents++;
    
    // trigger
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerSignal") && !passedHLT ) continue;   
    trigger++;

    // best tag-probe pair = mee closest to Z mass
    float minpull = 100000.;
    float mass    = 1000.;
    float okmass  = 1000.;

    for(int iele1=0; iele1<nEle; iele1++) {
      TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);

      for(int iele2=iele1+1; iele2<nEle; iele2++) {
        TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);
	
	eleTot++;

        if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[iele1]) || !m_selection->passCut("etaEleAcc",etaEle[iele2]) ) ) continue;
	eleEta++;

        if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc",electron1.Pt()) || !m_selection->passCut("ptEleAcc",electron2.Pt()) ) ) continue;
	elePt++;
	
        mass = (electron1+electron2).M();
        m_Zmass->Fill(mass);
        float pull=fabs(mass-91.1876);
	if(pull < minpull) {
	  okmass  = mass;
          minpull = pull;
          electrons[0] = iele1;
          electrons[1] = iele2;
        }
      }
    }

    if (okmass<999) twoele++;

    // start the tag & probe
    if( m_selection->passCut("meeWindow",okmass) ) {

      if ( okmass>110 || okmass<60 ) cout << "BACO!" << endl;
      invmass++;
      
      for(int iele=0; iele<2; ++iele) {
        int tag=electrons[0], probe=electrons[1];
        if(iele==1) {
          tag=electrons[1]; 
          probe=electrons[0];
        }

        TLorentzVector probeP4(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);
        TLorentzVector tagP4(pxEle[tag],    pyEle[tag],  pzEle[tag],  energyEle[tag]);

	// various about probe
        int charge = chargeEle[probe];
        float pt   = probeP4.Pt();
        float eta  = etaEle[probe];

        /// define the bins in which can be splitted the PDFs
        int iecal  = (fabs(etaEle[probe])<1.479) ? 0 : 1;
        int iptbin = (probeP4.Pt()<15.0) ? 0 : 1;
        int fullclassRaw = classificationEle[probe];
        int iclass     = -1;
        int ifullclass = -1;
        if      ( fullclassRaw == GOLDEN )    { iclass = 0; ifullclass = 0; }
        else if ( fullclassRaw == BIGBREM )   { iclass = 0; ifullclass = 1; }
        else if ( fullclassRaw == NARROW )    { iclass = 0; ifullclass = 2; }
        else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
        if (iclass>-1) tagProbeTot++;

        // apply the electron ID loose on the tag electron
	// float tagIdentified = anaUtils.electronIdVal(eleIdCutsEle[tag],eleIdRobustLoose);  
	bool tagIdentified, tagIsolated, tagConvRej;
        tagIdentified = tagIsolated = tagConvRej = false;
        isEleID(tag,&tagIdentified,&tagIsolated,&tagConvRej);
        if (tagIdentified && tagConvRej) tagId++;

        // apply tracker isolation on the tag electron 
	if (tagIsolated) tagIsol++;

        // apply tracker isolation on the probe electron
        bool probeIdentified, probeIsolated, probeConvRej;
        probeIdentified = probeIsolated = probeConvRej = true;
        if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
          isEleID(probe,&probeIdentified,&probeIsolated,&probeConvRej);
        }
	if(probeIsolated) probeIsol++;


	// some eleID variables
        int sc = superClusterIndexEle[probe];
        float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaSC[probe]));
        float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiSC[probe]));
        float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiSC[probe]));
        float s1s9          = eMaxSC[sc]/e3x3SC[sc];
        float s9s25         = e3x3SC[sc]/e5x5SC[sc];
	double dPhiCalo     = deltaPhiAtCaloEle[probe];
	double dPhiVtx      = deltaPhiAtVtxEle[probe];
	double dEtaVtx      = deltaEtaAtVtxEle[probe];
	double EoPout       = eSeedOverPoutEle[probe];
	double EoP          = eSuperClusterOverPEle[probe];
	double HoE          = hOverEEle[probe];
	int nbrem           = nbremsEle[probe];
        float fbrem         = fbremEle[probe];
        
	
        /// fill the electron ID pdfs only if:
        /// the tag is loose isolated and identified (ALWAYS)
        /// the probe is loose isolated              (ONLY IF REQUIRED)
        if( tagIsolated && tagIdentified && probeIsolated && iclass>-1) {   

          dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
          dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
          dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEtaVtx );
          EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
          HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
          sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
          sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
          sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
          s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
          s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );

          dPhiCaloClassEle      [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
          dPhiVtxClassEle       [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
          dEtaClassEle          [iecal][iptbin][iclass] -> Fill ( dEtaVtx );
          EoPoutClassEle        [iecal][iptbin][iclass] -> Fill ( EoPout );
          HoEClassEle           [iecal][iptbin][iclass] -> Fill ( HoE );
          sigmaIEtaIEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
          sigmaIEtaIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
          sigmaIPhiIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
          s1s9ClassEle          [iecal][iptbin][iclass] -> Fill ( s1s9 );
          s9s25ClassEle         [iecal][iptbin][iclass] -> Fill ( s9s25 );

          dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
          dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
          dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEtaVtx );
          EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
          HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
          sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
          sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
          sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
          s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
          s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );
	  
          // fill the reduced tree
          // changed: to be adapted
	  // reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,sigmaIEtaIEta,sigmaIPhiIPhi,fbrem,nbrem,pt,eta,charge);
          reducedTree.fillAttributesSignal(okmass,1,-1,-1,-1);
          reducedTree.fillCategories(iecal,iptbin,iclass,nbrem);
          reducedTree.store();

        } // fill histograms
	
      } // loop over the 2 Z electrons
      
    } // end tag and probe
    
  } // loop over events
  
  cout << "statistics from Tag and Probe: " << endl;
  cout << "allevents   = " << allevents << endl;
  cout << "trigger     = " << trigger << endl;
  cout << "twoele      = " << twoele  << endl;
  cout << "invmass     = " << invmass << endl;
  cout << "tagProbeTot = " << tagProbeTot << endl;
  cout << "tagId       = " << tagId       << endl;
  cout << "tagIsol     = " << tagIsol     << endl;
  cout << "probeIsol   = " << probeIsol   << endl;  
  cout << "statistics from Tag and Probe - electrons: " << endl;
  cout << "eleTot      = " << eleTot  << endl;
  cout << "eleEta      = " << eleEta  << endl;
  cout << "elePt       = " << elePt   << endl;

  reducedTree.save();
}


// PDF for the signal from a pure Z MC sample. Same requests as above
void LHPdfsProducer::LoopZ(const char *treefilesuffix) {

  if(fChain == 0) return;

  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesSignal();
  reducedTree.addCategories();
  
  // counters
  int allevents    = 0;
  int taus         = 0;
  int mc           = 0;
  int trigger      = 0;
  int twoele       = 0;
  int foundReco1   = 0;
  int foundReco2   = 0;
  int eleEta1      = 0;
  int elePt1       = 0;
  int tagProbeTot1 = 0;
  int probeIsol1   = 0;
  int eleEta2      = 0;
  int elePt2       = 0;
  int tagProbeTot2 = 0;
  int probeIsol2   = 0;

  // loop over events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();

  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    allevents++;

    // the stable particle list is truncated, if there is a tau not possible to say what happens... 
    bool tauPresence=false;
    for(int iMc=0; iMc<50; iMc++) {
      if ( (fabs(idMc[iMc])==15) ) { tauPresence=true; break; }
    }
    taus++;

    // to find the real electrons from Z
    int mcInd1 = -1;
    int mcInd2 = -1;
    for(int iMc=0; iMc<nMc; iMc++) {
      if ( (fabs(idMc[iMc])==11) && (fabs(idMc[mothMc[iMc]])==23) && mcInd1==-1 )              {mcInd1=iMc; continue;}
      if ( (fabs(idMc[iMc])==11) && (fabs(idMc[mothMc[iMc]])==23) && mcInd1!=-1 && mcInd2==-1) {mcInd2=iMc; break;}
    }
    TVector3 _mcEle1(0,0,0);
    TVector3 _mcEle2(0,0,0);
    if(mcInd1>-1) _mcEle1 = TVector3(pMc[mcInd1]*cos(phiMc[mcInd1])*sin(thetaMc[mcInd1]),pMc[mcInd1]*sin(phiMc[mcInd1])*sin(thetaMc[mcInd1]),pMc[mcInd1]*cos(thetaMc[mcInd1]));
    if(mcInd2>-1) _mcEle2 = TVector3(pMc[mcInd2]*cos(phiMc[mcInd2])*sin(thetaMc[mcInd2]),pMc[mcInd2]*sin(phiMc[mcInd2])*sin(thetaMc[mcInd2]),pMc[mcInd2]*cos(thetaMc[mcInd2]));
    if (mcInd1<0 || mcInd2<0) continue; 
    mc++;

    // this is for MC only. No json
    bool newTriggerMask = false;
    reloadTriggerMask(newTriggerMask);
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerSignal") && !passedHLT ) continue;   
    trigger++;

    // electrons matching MC truth
    float deltaRmin_mc1   = 999.;
    float deltaRmin_mc2   = 999.;
    int theClosestEle_mc1 = -1;
    int theClosestEle_mc2 = -1;
    for(int theEle=0; theEle<nEle; theEle++){
      
      TVector3 _p3Ele(pxEle[theEle],pyEle[theEle],pzEle[theEle]);
      float deltaR_mc1 = _p3Ele.DeltaR(_mcEle1);
      float deltaR_mc2 = _p3Ele.DeltaR(_mcEle2);
      
      if (deltaR_mc1<deltaRmin_mc1 && deltaR_mc1<0.5){
	deltaRmin_mc1     = deltaR_mc1;
	theClosestEle_mc1 = theEle;
      }

      if (deltaR_mc2<deltaRmin_mc2 && deltaR_mc2<0.5){
	deltaRmin_mc2     = deltaR_mc2;
	theClosestEle_mc2 = theEle;
      }
    }
    

    // fill PDFs using the electrons in the acceptance and matching the MC truth
    if (theClosestEle_mc1>-1) {
      
      foundReco1++;

      TLorentzVector electron(pxEle[theClosestEle_mc1],pyEle[theClosestEle_mc1],pzEle[theClosestEle_mc1],energyEle[theClosestEle_mc1]);
      if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[theClosestEle_mc1])) ) continue;
      eleEta1++;
      if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc", electron.Pt())) ) continue;
      elePt1++;

      // various
      int charge = chargeEle[theClosestEle_mc1];
      float pt   = electron.Pt();
      float eta  = etaEle[theClosestEle_mc1];

      // define the bins in which can be splitted the PDFs
      int iecal  = (fabs( etaEle[theClosestEle_mc1])<1.479) ? 0 : 1;
      int iptbin = (pt<15.0) ? 0 : 1;
      int fullclassRaw = classificationEle[theClosestEle_mc1];
      int iclass = -1;
      int ifullclass = -1;
      if      ( fullclassRaw == GOLDEN )    { iclass = 0; ifullclass = 0; }
      else if ( fullclassRaw == BIGBREM )   { iclass = 0; ifullclass = 1; }
      else if ( fullclassRaw == NARROW )    { iclass = 0; ifullclass = 2; }
      else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
      if (iclass>-1) tagProbeTot1++;

      // apply loose tracker isolation on the first electron
      bool isolated1 = true;
      float relativeIsol1 = dr04TkSumPtEle[theClosestEle_mc1]/pt;	
      if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
	isolated1 = ( m_selection->passCut("relSumPtTracks",relativeIsol1) );
      }
      probeIsol1++;

      // some eleID variables
      int sc = superClusterIndexEle[theClosestEle_mc1];
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaSC[sc]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiSC[sc]));
      float s1s9          = eMaxSC[sc]/e3x3SC[sc];
      float s9s25         = e3x3SC[sc]/e5x5SC[sc];
      double dPhiVtx      = deltaPhiAtVtxEle[theClosestEle_mc1];
      double dEtaVtx      = deltaEtaAtVtxEle[theClosestEle_mc1];
      double EoPout       = eSeedOverPoutEle[theClosestEle_mc1];
      double EoP          = eSuperClusterOverPEle[theClosestEle_mc1];
      double HoE          = hOverEEle[theClosestEle_mc1];
      int nbrem           = nbremsEle[theClosestEle_mc1];
      float fbrem         = fbremEle[theClosestEle_mc1];

      // fill the reduced tree     
      if( isolated1 && iclass>-1 ) {
        //	reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,sigmaIEtaIEta,sigmaIPhiIPhi,fbrem,nbrem,pt,eta,charge);
	reducedTree.fillAttributesSignal(9999.,1,-1,-1,-1);
	reducedTree.fillCategories(iecal,iptbin,iclass,nbrem);
	reducedTree.store();
      } // fill tree
    } // ok 1st electron


    // fill PDFs using the electrons in the acceptance and matching the MC truth
    if (theClosestEle_mc2>-1) {
      
      foundReco2++;

      TLorentzVector electron(pxEle[theClosestEle_mc2],pyEle[theClosestEle_mc2],pzEle[theClosestEle_mc2],energyEle[theClosestEle_mc2]);
      if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[theClosestEle_mc2])) ) continue;
      eleEta2++;
      if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc", electron.Pt())) ) continue;
      elePt2++;

      // various
      int charge = chargeEle[theClosestEle_mc2];
      float pt   = electron.Pt();
      float eta  = etaEle[theClosestEle_mc2];

      // define the bins in which can be splitted the PDFs
      int iecal  = (fabs( etaEle[theClosestEle_mc2])<1.479) ? 0 : 1;
      int iptbin = (pt<15.0) ? 0 : 1;
      int fullclassRaw = classificationEle[theClosestEle_mc2];
      int iclass = -1;
      int ifullclass = -1;
      if      ( fullclassRaw == GOLDEN )    { iclass = 0; ifullclass = 0; }
      else if ( fullclassRaw == BIGBREM )   { iclass = 0; ifullclass = 1; }
      else if ( fullclassRaw == NARROW )    { iclass = 0; ifullclass = 2; }
      else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
      if (iclass>-1) tagProbeTot2++;

      // apply loose tracker isolation on the first electron
      bool isolated2 = true;
      float relativeIsol2 = dr04TkSumPtEle[theClosestEle_mc2]/pt;	
      if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
	isolated2 = ( m_selection->passCut("relSumPtTracks",relativeIsol2) );
      }
      probeIsol2++;

      // some eleID variables
      int sc = superClusterIndexEle[theClosestEle_mc2];
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaSC[sc]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiSC[sc]));
      float s1s9          = eMaxSC[sc]/e3x3SC[sc];
      float s9s25         = e3x3SC[sc]/e5x5SC[sc];
      double dPhiVtx      = deltaPhiAtVtxEle[theClosestEle_mc2];
      double dEtaVtx      = deltaEtaAtVtxEle[theClosestEle_mc2];
      double EoPout       = eSeedOverPoutEle[theClosestEle_mc2];
      double EoP          = eSuperClusterOverPEle[theClosestEle_mc2];
      double HoE          = hOverEEle[theClosestEle_mc2];
      int nbrem           = nbremsEle[theClosestEle_mc2];
      float fbrem         = fbremEle[theClosestEle_mc2];

      // fill the reduced tree     
      if( isolated2 && iclass>-1 ) {
        //	reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,sigmaIEtaIEta,sigmaIPhiIPhi,fbrem,nbrem,pt,eta,charge);
	reducedTree.fillAttributesSignal(9999.,1,-1,-1,-1);
	reducedTree.fillCategories(iecal,iptbin,iclass,nbrem);
	reducedTree.store();
      } // fill tree
    } // ok 2nd electron
    


  } // loop over events
  

  cout << "statistics from MC: " << endl;
  cout << "allevents    = " << allevents    << endl;
  cout << "taus         = " << taus         << endl;  
  cout << "mc           = " << mc           << endl;  
  cout << "trigger      = " << trigger      << endl;
  cout << "foundReco1   = " << foundReco1   << endl;
  cout << "tagProbeTot1 = " << tagProbeTot1 << endl;
  cout << "probeIsol1   = " << probeIsol1   << endl;  
  cout << "foundReco2   = " << foundReco2   << endl;
  cout << "tagProbeTot2 = " << tagProbeTot2 << endl;
  cout << "probeIsol2   = " << probeIsol2   << endl;  

  cout << "statistics from MC - electrons: " << endl;
  cout << "eleEta1      = " << eleEta1      << endl;
  cout << "elePt1       = " << elePt1       << endl;
  cout << "eleEta2      = " << eleEta2      << endl;
  cout << "elePt2       = " << elePt2       << endl;

  reducedTree.save();
}

// PDF for electrons matching the MC truth from Z. To decouple from tag and probe we 
// apply the same T&P as for data driven method.
// Probe electrons are within the acceptance (and loose isolated in the tracker)
// The tag electron is the one with the best match to the Z
// the tag must be within the acceptance, tracker isolated and loose identified
void LHPdfsProducer::LoopZTagAndProbeForMcTruth(const char *treefilesuffix) {
  
  if(fChain == 0) return;

  bookHistos();

  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesSignal();
  reducedTree.addCategories();
  
  // counters
  int allevents   = 0;
  int trigger     = 0;
  int twoele      = 0;
  int eleTot      = 0;
  int eleEta      = 0;
  int elePt       = 0;
  int invmass     = 0;
  int tagId       = 0;
  int tagIsol     = 0;
  int probeIsol   = 0;
  int tagProbeTot = 0;

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
    
    // to find the real electrons from Z (MC truth)
    int mcInd1 = -1;
    int mcInd2 = -1;
    for(int iMc=0; iMc<nMc; iMc++) {
      if ( (fabs(idMc[iMc])==11) && (fabs(idMc[mothMc[iMc]])==23) && mcInd1==-1 )              { mcInd1=iMc; continue;}
      if ( (fabs(idMc[iMc])==11) && (fabs(idMc[mothMc[iMc]])==23) && mcInd1!=-1 && mcInd2==-1) { mcInd2=iMc; break; }
    }
    TVector3 _mcEle1(0,0,0);
    TVector3 _mcEle2(0,0,0);
    if(mcInd1>-1) _mcEle1 = TVector3(pMc[mcInd1]*cos(phiMc[mcInd1])*sin(thetaMc[mcInd1]),pMc[mcInd1]*sin(phiMc[mcInd1])*sin(thetaMc[mcInd1]),pMc[mcInd1]*cos(thetaMc[mcInd1]));
    if(mcInd2>-1) _mcEle2 = TVector3(pMc[mcInd2]*cos(phiMc[mcInd2])*sin(thetaMc[mcInd2]),pMc[mcInd2]*sin(phiMc[mcInd2])*sin(thetaMc[mcInd2]),pMc[mcInd2]*cos(thetaMc[mcInd2]));

    // this is for MC only: no json
    bool newTriggerMask = false;
    reloadTriggerMask(newTriggerMask);
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerSignal") && !passedHLT ) continue;   
    trigger++;
    
    // best tag-probe pair = mee closest to Z mass
    float minpull = 100000.;
    float mass    = 1000.;
    float okmass  = 1000.;

    for(int iele1=0; iele1<nEle; iele1++) {
      TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);

      for(int iele2=iele1+1; iele2<nEle; iele2++) {
        TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);
	
	eleTot++;

        if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[iele1]) || !m_selection->passCut("etaEleAcc",etaEle[iele2]) ) ) continue;
	eleEta++;

        if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc",electron1.Pt()) || !m_selection->passCut("ptEleAcc",electron2.Pt()) ) ) continue;
	elePt++;
	
        mass = (electron1+electron2).M();
        float pull=fabs(mass-91.1876);
	if(pull < minpull) {
	  okmass  = mass;
          minpull = pull;
          electrons[0] = iele1;
          electrons[1] = iele2;
        }
      }
    }

    if (okmass<999) twoele++;

    // start the tag & probe
    if( m_selection->passCut("meeWindow",okmass) ) {

      if ( okmass>110 || okmass<60 ) cout << "BACO!" << endl;
      invmass++;
      
      for(int iele=0; iele<2; ++iele) {
        int tag=electrons[0], probe=electrons[1];
        if(iele=1) {
          tag=electrons[1]; 
          probe=electrons[0];
        }
	
        TLorentzVector probeP4(pxEle[probe],pyEle[probe],pzEle[probe],energyEle[probe]);
        TLorentzVector tagP4(pxEle[tag],    pyEle[tag],  pzEle[tag],  energyEle[tag]);

	// various about probe
        int charge = chargeEle[probe];
        float pt   = probeP4.Pt();
        float eta  = etaEle[probe];

        /// define the bins in which can be splitted the PDFs
        int iecal  = (fabs(etaEle[probe])<1.479) ? 0 : 1;
        int iptbin = (probeP4.Pt()<15.0) ? 0 : 1;
        int fullclassRaw = classificationEle[probe];
        int iclass     = -1;
        int ifullclass = -1;
        if      ( fullclassRaw == GOLDEN )    { iclass = 0; ifullclass = 0; }
        else if ( fullclassRaw == BIGBREM )   { iclass = 0; ifullclass = 1; }
        else if ( fullclassRaw == NARROW )    { iclass = 0; ifullclass = 2; }
        else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
        if (iclass>-1) tagProbeTot++;

        // apply the electron ID loose on the tag electron
	// float tagIdentified = anaUtils.electronIdVal(eleIdCutsEle[tag],eleIdRobustLoose);  
        bool tagIdentified, tagIsolated, tagConvRej;
        tagIdentified = tagIsolated = tagConvRej = false;
        isEleID(tag,&tagIdentified,&tagIsolated,&tagConvRej);
        if (tagIdentified && tagConvRej) tagId++;

        // apply tracker isolation on the tag electron 
	if (tagIsolated) tagIsol++;

        // apply tracker isolation on the probe electron
        bool probeIdentified, probeIsolated, probeConvRej;
        probeIdentified = probeIsolated = probeConvRej = true;
        if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
          isEleID(probe,&probeIdentified,&probeIsolated,&probeConvRej);
        }
	if(probeIsolated) probeIsol++;


	// some eleID variables
        int sc = superClusterIndexEle[probe];
        float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaSC[sc]));
        float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiSC[sc]));
        float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiSC[sc]));
        float s1s9          = eMaxSC[sc]/e3x3SC[sc];
        float s9s25         = e3x3SC[sc]/e5x5SC[sc];
	double dPhiCalo     = deltaPhiAtCaloEle[probe];
	double dPhiVtx      = deltaPhiAtVtxEle[probe];
	double dEtaVtx      = deltaEtaAtVtxEle[probe];
	double EoPout       = eSeedOverPoutEle[probe];
	double EoP          = eSuperClusterOverPEle[probe];
	double HoE          = hOverEEle[probe];
	int nbrem           = nbremsEle[probe];
        float fbrem         = fbremEle[probe];

        /// fill the electron ID pdfs only if:
        /// the tag is loose isolated and identified (ALWAYS)
        /// the probe is loose isolated              (ONLY IF REQUIRED)
	/// it matches one of the MC electrons from Z
	if( tagIsolated && tagIdentified && probeIsolated && iclass>-1) {   

	  // does the probe match the MC truth?
	  bool matchMc = false;
	  TVector3 p3Probe(pxEle[probe],pyEle[probe],pzEle[probe]);
	  if (mcInd1>-1) { float deltaR = p3Probe.DeltaR(_mcEle1); if (deltaR<0.3) matchMc = true; }
	  if (mcInd2>-1) { float deltaR = p3Probe.DeltaR(_mcEle2); if (deltaR<0.3) matchMc = true; }

	  // only if it matches with MC fill the reduced tree
	  if (matchMc) { 
            //	    reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,sigmaIEtaIEta,sigmaIPhiIPhi,fbrem,nbrem,pt,eta,charge);
	    reducedTree.fillAttributesSignal(okmass,1,-1,-1,-1);
	    reducedTree.fillCategories(iecal,iptbin,iclass,nbrem);
	    reducedTree.store();

            dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
            dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
            dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEtaVtx );
            EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
            HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
            sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
            sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
            sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
            s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
            s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );

            dPhiCaloClassEle      [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
            dPhiVtxClassEle       [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
            dEtaClassEle          [iecal][iptbin][iclass] -> Fill ( dEtaVtx );
            EoPoutClassEle        [iecal][iptbin][iclass] -> Fill ( EoPout );
            HoEClassEle           [iecal][iptbin][iclass] -> Fill ( HoE );
            sigmaIEtaIEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
            sigmaIEtaIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
            sigmaIPhiIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
            s1s9ClassEle          [iecal][iptbin][iclass] -> Fill ( s1s9 );
            s9s25ClassEle         [iecal][iptbin][iclass] -> Fill ( s9s25 );

            dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
            dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
            dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEtaVtx );
            EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
            HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
            sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
            sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
            sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
            s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
            s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );

	  }

        } // fill histograms
	
      } // loop over the 2 Z electrons
      
    } // end tag and probe
    
  } // loop over events
  
  cout << "statistics from Tag and Probe: " << endl;
  cout << "allevents   = " << allevents << endl;
  cout << "trigger     = " << trigger << endl;
  cout << "twoele      = " << twoele  << endl;
  cout << "invmass     = " << invmass << endl;
  cout << "tagProbeTot = " << tagProbeTot << endl;
  cout << "tagId       = " << tagId       << endl;
  cout << "tagIsol     = " << tagIsol     << endl;
  cout << "probeIsol   = " << probeIsol   << endl;  
  cout << "statistics from Tag and Probe - electrons: " << endl;
  cout << "eleTot      = " << eleTot  << endl;
  cout << "eleEta      = " << eleEta  << endl;
  cout << "elePt       = " << elePt   << endl;

  reducedTree.save();
}


// PDF for the signal from a pure Z MC sample. Same requests as above
void LHPdfsProducer::LoopZwithMass(const char *treefilesuffix) {

  if(fChain == 0) return;

  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesSignal();
  reducedTree.addCategories();
  
  // counters
  int allevents    = 0;
  int mc           = 0;
  int trigger      = 0;
  int twoele       = 0;
  int foundReco    = 0;
  int eleEta       = 0;
  int elePt        = 0;
  int tagProbeTot  = 0;
  int probeIsol    = 0;

  // loop over events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();

  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    allevents++;

    // to find the real electrons from Z
    int mcInd1 = -1;
    int mcInd2 = -1;
    for(int iMc=0; iMc<nMc; iMc++) {
      if ( (fabs(idMc[iMc])==11) && (fabs(idMc[mothMc[iMc]])==23) && mcInd1==-1 )              { mcInd1=iMc; continue;}
      if ( (fabs(idMc[iMc])==11) && (fabs(idMc[mothMc[iMc]])==23) && mcInd1!=-1 && mcInd2==-1) { mcInd2=iMc; break; }
    }
    TVector3 _mcEle1(0,0,0);
    TVector3 _mcEle2(0,0,0);
    if(mcInd1>-1) _mcEle1 = TVector3(pMc[mcInd1]*cos(phiMc[mcInd1])*sin(thetaMc[mcInd1]),pMc[mcInd1]*sin(phiMc[mcInd1])*sin(thetaMc[mcInd1]),pMc[mcInd1]*cos(thetaMc[mcInd1]));
    if(mcInd2>-1) _mcEle2 = TVector3(pMc[mcInd2]*cos(phiMc[mcInd2])*sin(thetaMc[mcInd2]),pMc[mcInd2]*sin(phiMc[mcInd2])*sin(thetaMc[mcInd2]),pMc[mcInd2]*cos(thetaMc[mcInd2]));
    if (mcInd1<0 || mcInd2<0) continue; 
    mc++;

    // this is for MC only: no json
    bool newTriggerMask = false;
    reloadTriggerMask(newTriggerMask);
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerSignal") && !passedHLT ) continue;   
    trigger++;

    // electrons matching MC truth
    float deltaRmin_mc1   = 999.;
    float deltaRmin_mc2   = 999.;
    int theClosestEle_mc1 = -1;
    int theClosestEle_mc2 = -1;
    for(int theEle=0; theEle<nEle; theEle++){
      
      TVector3 _p3Ele(pxEle[theEle],pyEle[theEle],pzEle[theEle]);
      float deltaR_mc1 = _p3Ele.DeltaR(_mcEle1);
      float deltaR_mc2 = _p3Ele.DeltaR(_mcEle2);
      
      if (deltaR_mc1<deltaRmin_mc1 && deltaR_mc1<0.5){
	deltaRmin_mc1     = deltaR_mc1;
	theClosestEle_mc1 = theEle;
      }

      if (deltaR_mc2<deltaRmin_mc2 && deltaR_mc2<0.5){
	deltaRmin_mc2     = deltaR_mc2;
	theClosestEle_mc2 = theEle;
      }
    }
    

    // fill PDFs using the electrons in the acceptance and matching the MC truth
    if (theClosestEle_mc1>-1 && theClosestEle_mc2>-1) {
      
      foundReco++;

      TLorentzVector electron1(pxEle[theClosestEle_mc1],pyEle[theClosestEle_mc1],pzEle[theClosestEle_mc1],energyEle[theClosestEle_mc1]);
      TLorentzVector electron2(pxEle[theClosestEle_mc2],pyEle[theClosestEle_mc2],pzEle[theClosestEle_mc2],energyEle[theClosestEle_mc2]);
      if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[theClosestEle_mc1])) ) continue;
      if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[theClosestEle_mc2])) ) continue;
      eleEta++;
      if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc", electron1.Pt())) ) continue;
      if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc", electron2.Pt())) ) continue;
      elePt++;

      // various
      int charge1 = chargeEle[theClosestEle_mc1];
      float pt1   = electron1.Pt();
      float eta1  = etaEle[theClosestEle_mc1];
      int charge2 = chargeEle[theClosestEle_mc2];
      float pt2   = electron2.Pt();
      float eta2  = etaEle[theClosestEle_mc2];

      // define the bins in which can be splitted the PDFs
      int iecal1  = (fabs( etaEle[theClosestEle_mc1])<1.479) ? 0 : 1;
      int iecal2  = (fabs( etaEle[theClosestEle_mc2])<1.479) ? 0 : 1;
      int iptbin1 = (pt1<15.0) ? 0 : 1;
      int iptbin2 = (pt2<15.0) ? 0 : 1;
      int fullclassRaw1 = classificationEle[theClosestEle_mc1];
      int fullclassRaw2 = classificationEle[theClosestEle_mc2];
      int iclass1 = -1;
      int ifullclass1 = -1;
      int iclass2 = -1;
      int ifullclass2 = -1;
      if      ( fullclassRaw1 == GOLDEN )    { iclass1 = 0; ifullclass1 = 0; }
      else if ( fullclassRaw1 == BIGBREM )   { iclass1 = 0; ifullclass1 = 1; }
      else if ( fullclassRaw1 == NARROW )    { iclass1 = 0; ifullclass1 = 2; }
      else if ( fullclassRaw1 == SHOWERING ) { iclass1 = 1; ifullclass1 = 3; }
      if      ( fullclassRaw2 == GOLDEN )    { iclass2 = 0; ifullclass2 = 0; }
      else if ( fullclassRaw2 == BIGBREM )   { iclass2 = 0; ifullclass2 = 1; }
      else if ( fullclassRaw2 == NARROW )    { iclass2 = 0; ifullclass2 = 2; }
      else if ( fullclassRaw2 == SHOWERING ) { iclass2 = 1; ifullclass2 = 3; }
      if (iclass1>-1 && iclass2>-1) tagProbeTot++;

      // apply loose tracker isolation 
      bool isolated1 = true;
      bool isolated2 = true;
      float relativeIsol1 = dr04TkSumPtEle[theClosestEle_mc1]/pt1;	
      float relativeIsol2 = dr04TkSumPtEle[theClosestEle_mc2]/pt2;	
      if ( m_selection->getSwitch("applyIsolationOnProbe") ) {
	isolated1 = ( m_selection->passCut("relSumPtTracks",relativeIsol1) );
	isolated2 = ( m_selection->passCut("relSumPtTracks",relativeIsol2) );
      }
      if (isolated1 && isolated2) probeIsol++;

      // some eleID variables
      int sc1 = superClusterIndexEle[theClosestEle_mc1];
      float sigmaIEtaIEta1 = sqrt(fabs(covIEtaIEtaSC[sc1]));
      float sigmaIPhiIPhi1 = sqrt(fabs(covIPhiIPhiSC[sc1]));
      float s1s91          = eMaxSC[sc1]/e3x3SC[sc1];
      float s9s251         = e3x3SC[sc1]/e5x5SC[sc1];
      double dPhiVtx1      = deltaPhiAtVtxEle[theClosestEle_mc1];
      double dEtaVtx1      = deltaEtaAtVtxEle[theClosestEle_mc1];
      double EoPout1       = eSeedOverPoutEle[theClosestEle_mc1];
      double EoP1          = eSuperClusterOverPEle[theClosestEle_mc1];
      double HoE1          = hOverEEle[theClosestEle_mc1];
      int nbrem1           = nbremsEle[theClosestEle_mc1];
      float fbrem1         = fbremEle[theClosestEle_mc1];
      int sc2 = superClusterIndexEle[theClosestEle_mc2];
      float sigmaIEtaIEta2 = sqrt(fabs(covIEtaIEtaSC[sc2]));
      float sigmaIPhiIPhi2 = sqrt(fabs(covIPhiIPhiSC[sc2]));
      float s1s92          = eMaxSC[sc2]/e3x3SC[sc2];
      float s9s252         = e3x3SC[sc2]/e5x5SC[sc2];
      double dPhiVtx2      = deltaPhiAtVtxEle[theClosestEle_mc2];
      double dEtaVtx2      = deltaEtaAtVtxEle[theClosestEle_mc2];
      double EoPout2       = eSeedOverPoutEle[theClosestEle_mc2];
      double EoP2          = eSuperClusterOverPEle[theClosestEle_mc2];
      double HoE2          = hOverEEle[theClosestEle_mc2];
      int nbrem2           = nbremsEle[theClosestEle_mc2];
      float fbrem2         = fbremEle[theClosestEle_mc2];

      float mass = (electron1+electron2).M();
      if( m_selection->passCut("meeWindow",mass) ) {

	// fill the reduced tree     
	if( isolated1 && isolated2 && iclass1>-1 && iclass2>-1 ) {
          //	  reducedTree.fillVariables(EoPout1,EoP1,HoE1,dEtaVtx1,dPhiVtx1,s9s251,sigmaIEtaIEta1,sigmaIPhiIPhi1,fbrem1,nbrem1,pt1,eta1,charge1);
	  reducedTree.fillAttributesSignal(mass,1,-1,-1,-1);
	  reducedTree.fillCategories(iecal1,iptbin1,iclass1,nbrem1);
          //	  reducedTree.fillVariables(EoPout2,EoP2,HoE2,dEtaVtx2,dPhiVtx2,s9s252,sigmaIEtaIEta2,sigmaIPhiIPhi2,fbrem2,nbrem2,pt2,eta2,charge2);
	  reducedTree.fillAttributesSignal(mass,1,-1,-1,-1);
	  reducedTree.fillCategories(iecal2,iptbin2,iclass2,nbrem2);
	  
	  reducedTree.store();
	} // fill tree
      } // mass
    } // ok 2 electron
  } // loop over events
  

  cout << "statistics from MC: " << endl;
  cout << "allevents    = " << allevents   << endl;
  cout << "mc           = " << mc          << endl;  
  cout << "trigger      = " << trigger     << endl;
  cout << "foundReco    = " << foundReco   << endl;
  cout << "tagProbeTot  = " << tagProbeTot << endl;
  cout << "probeIsol    = " << probeIsol   << endl;  

  cout << "statistics from MC - electrons: " << endl;
  cout << "eleEta      = " << eleEta         << endl;
  cout << "elePt       = " << elePt          << endl;

  reducedTree.save();
}

void LHPdfsProducer::LoopQCD() {
  
  if(fChain == 0) return;
  
  bookHistos();

  // json 
  unsigned int lastLumi=0;
  unsigned int lastRun=0;
  
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // reload trigger mask
    bool newTriggerMask = false;
    if(isData_) newTriggerMask = true;
    reloadTriggerMask(newTriggerMask);
    
    // Good Run selection 
    if (isData_ && m_selection->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun  = runNumber;
        lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }    
    if (isData_ && m_selection->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }

    // trigger
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerQCDBack") && !passedHLT ) continue;   

    // fill the PDFs for QCD with all the (isolated) reco'ed electrons
    for(int iele=0;iele<nEle;iele++) {
      
      if( m_selection->getSwitch("etaEleAcc") && 
          ! m_selection->passCut("etaEleAcc",etaEle[iele]) ) continue;
      
      if ( m_selection->getSwitch("applyIsolationOnProbe") &&
           ! m_selection->passCut("relSumPtTracks",dr04TkSumPtEle[iele]) ) continue;
      
      TLorentzVector eleP4(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);

      /// define the bins in which can be splitted the PDFs
      int iecal  = (fabs( etaEle[iele])<1.479) ? 0 : 1;
      int iptbin = (eleP4.Pt()<15.0) ? 0 : 1;
      
      int fullclassRaw = classificationEle[iele];
        
      int iclass = -1;
      int ifullclass = -1;
      if ( fullclassRaw == 0 || fullclassRaw == 100 ) { // golden
        iclass = 0;
        ifullclass = 0;
      }
      else if ( fullclassRaw == 10 || fullclassRaw == 110 ) { // bigbrem
        iclass = 0;
        ifullclass = 1;
      }
      else if ( fullclassRaw == 20 || fullclassRaw == 120 ) { // narrow
        iclass = 0;
        ifullclass = 2;
      }
      else if ( (fullclassRaw >= 30 && fullclassRaw <= 40) ||
                (fullclassRaw >= 130 && fullclassRaw <= 140) ) { // showering + cracks
        iclass = 1;
        ifullclass = 3;
      }
        
      int sc = superClusterIndexEle[iele];
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaSC[sc]));
      float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiSC[sc]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiSC[sc]));
      float s1s9   = eMaxSC[sc]/e3x3SC[sc];
      float s9s25  = e3x3SC[sc]/e5x5SC[sc];
      double dPhiCalo = deltaPhiAtCaloEle[iele];
      double dPhiVtx  = deltaPhiAtVtxEle[iele];
      double dEtaVtx  = deltaEtaAtVtxEle[iele];
      double EoPout   = eSeedOverPoutEle[iele];
      double EoP      = eSuperClusterOverPEle[iele];
      double HoE = hOverEEle[iele];

      dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEtaVtx );
      EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
      sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
      s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );

      dPhiCaloClassEle        [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle         [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle            [iecal][iptbin][iclass] -> Fill ( dEtaVtx );
      EoPoutClassEle          [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle             [iecal][iptbin][iclass] -> Fill ( HoE );
      sigmaIEtaIEtaClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
      s1s9ClassEle            [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle           [iecal][iptbin][iclass] -> Fill ( s9s25 );

      dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEtaVtx );
      EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
      sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
      s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );

    } // loop on electrons

  } // loop on events
  
}


void LHPdfsProducer::LoopQCDTagAndProbe(const char *treefilesuffix) {

  if(fChain == 0) return;
  
  bookHistos();

  char treename[200];
  sprintf(treename,"%s",treefilesuffix);
  RedEleIDTree reducedTree(treename);
  reducedTree.addAttributesBackground();
  reducedTree.addCategories();
  reducedTree.addGamma();        // to find the best cut for anti-isolationA

  // counters
  int nocrack = 0;

  // json 
  unsigned int lastLumi = 0;
  unsigned int lastRun  = 0;

  // loop over events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   
    nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // reload trigger mask
    bool newTriggerMask = false;
    if(isData_) newTriggerMask = true;
    reloadTriggerMask(newTriggerMask);
    
    // Good Run selection 
    if (isData_ && m_selection->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun  = runNumber;
        lastLumi = lumiBlock;
	std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    if (isData_ && m_selection->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    // count all events
    m_counters->IncrVar("allEvents");

    // pT hat cut
    if( m_selection->getSwitch("ptHat") && (!m_selection->passCut("ptHat",genPtHat) ) ) continue;
    m_counters->IncrVar("pthat");

    // QCD trigger 
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerQCDBack") && !passedHLT ) continue;   
    m_counters->IncrVar("trigger");

    // electrons and jets within acceptance
    vector<int> probeCandidates;
    vector<int> tagCandidates;
    
    // selecting electrons in the acceptance and possibly loose isolated to work as probe
    // spikes removal applied here
    for(int iele=0; iele<nEle; iele++) {
      TLorentzVector electron(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);
      float relativeIsol = dr04TkSumPtEle[iele]/electron.Pt();
      if( m_selection->getSwitch("etaEleAcc") && (!m_selection->passCut("etaEleAcc",etaEle[iele]) ) ) continue;
      if( m_selection->getSwitch("ptEleAcc")  && (!m_selection->passCut("ptEleAcc",electron.Pt()) ) ) continue;      
      if( m_selection->getSwitch("applyIsolationOnProbe") && !m_selection->passCut("relSumPtTracks",relativeIsol) ) continue;
      // removing spikes
      int theSuperCluster = superClusterIndexEle[iele];
      float e4SwissCross  = e4SwissCrossSC[theSuperCluster];
      float e1            = eMaxSC[theSuperCluster];
      if(m_selection->getSwitch("spikeFraction") && !m_selection->passCut("spikeFraction", 1.0-e4SwissCross/e1) ) continue;
      probeCandidates.push_back(iele);
    }
    if (probeCandidates.size()>0) m_counters->IncrVar("oneele");

    // selecting jets in the acceptance to work as tag
    for (int ijet=0; ijet<nAK5PFPUcorrJet; ijet++){    
      TVector3 p3Jet(pxAK5PFPUcorrJet[ijet],pyAK5PFPUcorrJet[ijet],pzAK5PFPUcorrJet[ijet]);
      if ( m_selection->getSwitch("etaJetAcc") && !m_selection->passCut("etaJetAcc",etaAK5PFPUcorrJet[ijet]) ) continue;
      if ( m_selection->getSwitch("etJetAcc")  && !m_selection->passCut("etJetAcc",p3Jet.Pt()) )   continue;
      tagCandidates.push_back(ijet);
    }
    if (tagCandidates.size()>0) m_counters->IncrVar("onejet"); 


    // to select the highest Et jet in the event within the acceptance in eta (tag) matching a probe electron
    int theTag   = -999;
    int theProbe = -999;
    float tagEt  = -999.;
    
    for (int iJet=0; iJet<tagCandidates.size(); iJet++){ 

      TLorentzVector p4Jet;
      TVector3 p3Jet(pxAK5PFPUcorrJet[tagCandidates[iJet]],pyAK5PFPUcorrJet[tagCandidates[iJet]],pzAK5PFPUcorrJet[tagCandidates[iJet]]);
      p4Jet.SetXYZT (pxAK5PFPUcorrJet[tagCandidates[iJet]],pyAK5PFPUcorrJet[tagCandidates[iJet]],pzAK5PFPUcorrJet[tagCandidates[iJet]], energyAK5PFPUcorrJet[tagCandidates[iJet]]);

      // we use the highest ET jet (with a potential probe) as a tag 
      if (p3Jet.Pt()<tagEt) continue;

      // we choose the probe with max Dphi wrt the tag
      float dPhiMax = -999.;

      for (int iEle=0; iEle<probeCandidates.size(); iEle++){       

	m_counters->IncrVar("eletot"); 

	TLorentzVector p4Ele;	
	TVector3 p3Ele(pxEle[probeCandidates[iEle]],pyEle[probeCandidates[iEle]],pzEle[probeCandidates[iEle]]);
	p4Ele.SetXYZT (pxEle[probeCandidates[iEle]],pyEle[probeCandidates[iEle]],pzEle[probeCandidates[iEle]],energyEle[probeCandidates[iEle]]);

	// the electron and the tag must not be the same thing
	float deltaR = p3Ele.DeltaR(p3Jet);
	if (deltaR<0.3) continue;

	// minimal requirements on invariant mass and separation
	float deltaPhi = fabs(p3Jet.DeltaPhi(p3Ele));
	float invMass  = (p4Jet+p4Ele).M();
	if ( m_selection->getSwitch("jetDeltaPhi") && !m_selection->passCut("jetDeltaPhi", deltaPhi) ) continue;
	m_counters->IncrVar("deltaphi"); 

	// we have a potential probe and this is highest ET jet up to now -> it's our tag
	tagEt  = p3Jet.Pt();
	theTag = tagCandidates[iJet];
	
	// in case of several probes for this tag
	if (deltaPhi>dPhiMax) { 
	  dPhiMax  = deltaPhi;
	  theProbe = probeCandidates[iEle];
	}
      }
    }
  
    // we need a tag and a probe
    if (theTag<-800 || theProbe<-800) continue;
    m_counters->IncrVar("tagandprobe"); 
      
    // variables for the tree
    TLorentzVector p4Tag, p4Probe;	
    TVector3 p3Probe(pxEle[theProbe],pyEle[theProbe],pzEle[theProbe]);
    p4Probe.SetXYZT (pxEle[theProbe],pyEle[theProbe],pzEle[theProbe],energyEle[theProbe]);
    TVector3 p3Tag(pxAK5PFPUcorrJet[theTag],pyAK5PFPUcorrJet[theTag],pzAK5PFPUcorrJet[theTag]);
    p4Tag.SetXYZT (pxAK5PFPUcorrJet[theTag],pyAK5PFPUcorrJet[theTag],pzAK5PFPUcorrJet[theTag], energyAK5PFPUcorrJet[theTag]);
    float theDeltaPhi = fabs(p3Tag.DeltaPhi(p3Probe));
    float theInvMass  = (p4Tag+p4Probe).M();
    TVector3 p3Met(pxMet[0],pyMet[0],0.0);
    float theMet      = p3Met.Pt();
    if ( m_selection->getSwitch("met")  && !m_selection->passCut("met", theMet) ) continue;
    m_counters->IncrVar("met"); 
    
    // reject events with >=2 electrons making a Z
    if ( nEle>=2 ) {
      float mass = -1;
      float meepull = 1000.;      
      for (int iele=0; iele<nEle && iele != theProbe; iele++) {
        TLorentzVector p4Electron(pxEle[iele],pyEle[iele],pzEle[iele],energyEle[iele]);
        float mee = (p4Probe + p4Electron).M();
        if ( fabs( mee - 91.1876) < meepull ) {
          meepull = fabs(mee - 91.1876);
          mass = mee;
        }
      }
      if ( m_selection->getSwitch("jetInvMass")  && m_selection->passCut("jetInvMass", mass) ) continue;
    }
    m_counters->IncrVar("invmass"); 

    // only for Gamma+jets: check if the probe is isolated and if it close to the MC gamma
    // to be commented out if not running on gamma+jets  
    int isGamma = 0;
    /* 
       int mcGamma = -1;
       for(int iMc=0; iMc<10; iMc++) {   
       if ( idMc[iMc]==22 && mcGamma<0) { mcGamma=iMc; }
       }
       if (mcGamma<0) { cout << "gamma not found" << endl; continue; }
       TVector3 p3McGamma(pMc[mcGamma]*cos(phiMc[mcGamma])*sin(thetaMc[mcGamma]),pMc[mcGamma]*sin(phiMc[mcGamma])*sin(thetaMc[mcGamma]),pMc[mcGamma]*cos(thetaMc[mcGamma]));
       float deltaR_GammaProbe = p3Probe.DeltaR(p3McGamma);
       if (deltaR_GammaProbe<0.3) isGamma = 1;
    */
    float absTrackerIsolGamma = dr04TkSumPtEle[theProbe];
    float absEcalIsolGamma    = dr04EcalRecHitSumEtEle[theProbe];  
    float absHcalIsolGamma    = dr04HcalTowerSumEtEle[theProbe];   

    // others
    int charge           = chargeEle[theProbe];
    float pt             = p4Probe.Pt();
    float eta            = etaEle[theProbe];
    float absTrackerIsol = dr04TkSumPtEle[theProbe];
    float absEcalIsol    = dr04EcalRecHitSumEtEle[theProbe];  

    // to reject gamma+jets events: looking for not isolated probes
    if ( m_selection->getSwitch("antiIsolTracker") && !m_selection->passCut("antiIsolTracker", absTrackerIsol) ) continue;
    m_counters->IncrVar("trackerNotIsol"); 
    if ( m_selection->getSwitch("antiIsolEcal") && !m_selection->passCut("antiIsolEcal", absEcalIsol) ) continue;
    m_counters->IncrVar("ecalNotIsol"); 

    // full selection
    m_counters->IncrVar("fullSelection");

    // fill the PDFs for QCD with the probe
    TLorentzVector eleP4(pxEle[theProbe],pyEle[theProbe],pzEle[theProbe],energyEle[theProbe]);

    /// define the bins in which can be splitted the PDFs
    int iecal  = (fabs( etaEle[theProbe])<1.479) ? 0 : 1;
    int iptbin = (eleP4.Pt()<15.0) ? 0 : 1;
    int fullclassRaw = classificationEle[theProbe];
    int iclass = -1;
    int ifullclass = -1;    
    if      ( fullclassRaw == GOLDEN ) { iclass = 0; ifullclass = 0; }
    else if ( fullclassRaw == BIGBREM ) { iclass = 0; ifullclass = 1; }
    else if ( fullclassRaw == NARROW ) { iclass = 0; ifullclass = 2; }
    else if ( fullclassRaw == SHOWERING ) { iclass = 1; ifullclass = 3; }
    if (iclass>-1) nocrack++;
          
    if (iclass>-1) {
      int sc = superClusterIndexEle[theProbe];
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaSC[sc]));
      float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiSC[sc]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiSC[sc]));
      float s1s9        = eMaxSC[sc]/e3x3SC[sc];
      float s9s25       = e3x3SC[sc]/e5x5SC[sc];
      double dPhiCalo   = deltaPhiAtCaloEle[theProbe];
      double dPhiVtx    = deltaPhiAtVtxEle[theProbe];
      double dEtaVtx    = deltaEtaAtVtxEle[theProbe];
      double EoPout     = eSeedOverPoutEle[theProbe];
      double EoP        = eSuperClusterOverPEle[theProbe];
      double HoE        = hOverEEle[theProbe];
      float fbrem       = fbremEle[theProbe];
      int nbrem         = nbremsEle[theProbe];
      
      // conversion rejection variables
      int theMatchedTrack   = trackIndexEle[theProbe];
      int theExpInnerLayers = expInnerLayersTrack[theMatchedTrack];
      float convDcot        = convDcotEle[theProbe];
      float convDist        = convDistEle[theProbe];
            
      // histos
      dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEtaVtx );
      EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
      sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
      s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );
      
      dPhiCaloClassEle        [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle         [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle            [iecal][iptbin][iclass] -> Fill ( dEtaVtx );
      EoPoutClassEle          [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle             [iecal][iptbin][iclass] -> Fill ( HoE );
      sigmaIEtaIEtaClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiClassEle   [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
      s1s9ClassEle            [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle           [iecal][iptbin][iclass] -> Fill ( s9s25 );
      
      dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEtaVtx );
      EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
      sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
      s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );

      // fill the reduced tree
      //      reducedTree.fillVariables(EoPout,EoP,HoE,dEtaVtx,dPhiVtx,s9s25,s1s9,sigmaIEtaIEta,sigmaIPhiIPhi,fbrem,theExpInnerLayers,convDcot,convDist,pt,eta,charge); 

      reducedTree.fillAttributesBackground(theDeltaPhi,theInvMass,theMet,genPtHat);

      reducedTree.fillCategories(iecal,iptbin,iclass,nbrem);

      reducedTree.fillGamma(absTrackerIsolGamma,absEcalIsolGamma,absHcalIsolGamma,isGamma);

      reducedTree.store();
    } // no crack

  } // loop on events

  reducedTree.save();    
}


// make jet PDFs from W+jets 
void LHPdfsProducer::LoopWjets() {

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

    // this is for MC only
    bool newTriggerMask = false;
    reloadTriggerMask(newTriggerMask);
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerQCDBack") && !passedHLT ) continue;   
    
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
      
      if ( m_selection->getSwitch("applyIsolationOnProbe") &&
           ! m_selection->passCut("relSumPtTracks",dr04TkSumPtEle[iele]) ) continue;
      
      
      // define the bins in which can be splitted the PDFs
      int iecal = (fabs( etaEle[iele])<1.479) ? 0 : 1;
      int iptbin = (eleP4.Pt()<15.0) ? 0 : 1;
      
      int fullclassRaw = classificationEle[iele];
        
      int iclass = -1;
      int ifullclass = -1;
      if ( fullclassRaw == 0 || fullclassRaw == 100 ) { // golden
        iclass = 0;
        ifullclass = 0;
      }
      else if ( fullclassRaw == 10 || fullclassRaw == 110 ) { // bigbrem
        iclass = 0;
        ifullclass = 1;
      }
      else if ( fullclassRaw == 20 || fullclassRaw == 120 ) { // narrow
        iclass = 0;
        ifullclass = 2;
      }
      else if ( (fullclassRaw >= 30 && fullclassRaw <= 40) ||
                (fullclassRaw >= 130 && fullclassRaw <= 140) ) { // showering + cracks
        iclass = 1;
        ifullclass = 3;
      }

      int sc = superClusterIndexEle[iele];
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaSC[sc]));
      float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiSC[sc]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiSC[sc]));
      float s1s9   = eMaxSC[sc]/e3x3SC[sc];
      float s9s25  = e3x3SC[sc]/e5x5SC[sc];
      double dPhiCalo = deltaPhiAtCaloEle[iele];
      double dPhiVtx  = deltaPhiAtVtxEle[iele];
      double dEta     = deltaEtaAtVtxEle[iele];
      double EoPout   = eSeedOverPoutEle[iele];
      double EoP      = eSuperClusterOverPEle[iele];
      double HoE      = hOverEEle[iele];
      //      double dxy      = eleTrackDxyEle[iele];
      //      double dxySig   = eleTrackDxyEle[iele]/eleTrackDxyErrorEle[iele];
      
      dPhiCaloUnsplitEle      [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle       [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle          [iecal][iptbin] -> Fill ( dEta );
      EoPoutUnsplitEle        [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle           [iecal][iptbin] -> Fill ( HoE );
      sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
      s1s9UnsplitEle          [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle         [iecal][iptbin] -> Fill ( s9s25 );

      dPhiCaloClassEle      [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle       [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle          [iecal][iptbin][iclass] -> Fill ( dEta );
      EoPoutClassEle        [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle           [iecal][iptbin][iclass] -> Fill ( HoE );
      sigmaIEtaIEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
      s1s9ClassEle          [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle         [iecal][iptbin][iclass] -> Fill ( s9s25 );

      dPhiCaloFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle       [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle          [iecal][iptbin][ifullclass] -> Fill ( dEta );
      EoPoutFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle           [iecal][iptbin][ifullclass] -> Fill ( HoE );
      sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
      s1s9FullclassEle          [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle         [iecal][iptbin][ifullclass] -> Fill ( s9s25 );

    } // loop on electrons

  } // loop on events
  
}


// make jets PDF removing the 2 electrons from Z: TOCHECK!!! something strange!
void LHPdfsProducer::LoopZjets(const char *outname) {

  int nbinsEta = 60;
  float minEta = -2.5;
  float maxEta = 2.5;
  
  TH1F *AllFakesEta = new TH1F("AllFakesEta", "all reconstructed fakes", nbinsEta, minEta, maxEta);
  TH1F *UnmatchedFakesEta = new TH1F("UnmatchedFakesEta", "Fake electrons unmatched with true electrons", nbinsEta, minEta, maxEta);

  int nbinsPt = 20;
  float minPt = 10.0;
  float maxPt = 100.;

  TH1F *AllFakesPt = new TH1F("AllFakesPt", "all reconstructed fakes", nbinsPt, minPt, maxPt);
  TH1F *UnmatchedFakesPt = new TH1F("UnmatchedFakesPt", "Fake electrons unmatched with true electrons", nbinsPt, minPt, maxPt);

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
    
    // this is pure MC...
    bool newTriggerMask = false;
    reloadTriggerMask(newTriggerMask);
    Utils anaUtils;
    bool passedHLT = hasPassedHLT();
    if ( m_selection->getSwitch("requireTriggerSignal") && !passedHLT ) continue;   

    bool tauPresence=false;
    bool spuriousElectronsMadgraph=false;
    int mcEleMinusIndex=-1, mcElePlusIndex=-1;
    for(int imc=0; imc<50; imc++) {
      // not only ele from Zee: there is enu emission (V_ud?) in madgraph: remove these spurious electrons
      if ( (fabs(idMc[imc])==11) && fabs(idMc[mothMc[imc]])!=23 && fabs(idMc[mothMc[imc]])!=11 ) {
        spuriousElectronsMadgraph=true;
        break;
      }
      if ( idMc[imc]==11 && fabs(idMc[mothMc[imc]])==23 ) mcEleMinusIndex = imc;
      if ( idMc[imc]==-11 && fabs(idMc[mothMc[imc]])==23 ) mcElePlusIndex = imc;
      // since the stable particle list is truncated, if there is a tau 
      // not possible to say what happens...
      if ( (fabs(idMc[imc])==15) ) {
        tauPresence=true;
        break;
      }
    }

    if(tauPresence || spuriousElectronsMadgraph) continue;

    TVector3 mcElePlus(0,0,0), mcEleMinus(0,0,0);
    if(mcEleMinusIndex>-1) mcEleMinus = TVector3(pMc[mcEleMinusIndex]*cos(phiMc[mcEleMinusIndex])*sin(thetaMc[mcEleMinusIndex]), 
                                                 pMc[mcEleMinusIndex]*sin(phiMc[mcEleMinusIndex])*sin(thetaMc[mcEleMinusIndex]),
                                                 pMc[mcEleMinusIndex]*cos(thetaMc[mcEleMinusIndex]));
    
    if(mcElePlusIndex>-1) mcEleMinus = TVector3(pMc[mcElePlusIndex]*cos(phiMc[mcElePlusIndex])*sin(thetaMc[mcElePlusIndex]), 
                                                pMc[mcElePlusIndex]*sin(phiMc[mcElePlusIndex])*sin(thetaMc[mcElePlusIndex]),
                                                pMc[mcElePlusIndex]*cos(thetaMc[mcElePlusIndex]));
    

    // fill the PDFs for jets excluding the real electron of Z->ee
    // the two electrons are the nes that give the |mee-mZ| closest to 0
    float minpullZee = 1000.;
    electrons[0]=-1; electrons[1]=-1;
    if(nEle>=2) {
      for(int iele1=0; iele1<nEle; iele1++) {
        TLorentzVector electron1(pxEle[iele1],pyEle[iele1],pzEle[iele1],energyEle[iele1]);
        for(int iele2=iele1+1; iele2<nEle; iele2++) {
          TLorentzVector electron2(pxEle[iele2],pyEle[iele2],pzEle[iele2],energyEle[iele2]);
          
          float mass = (electron1+electron2).M();
          m_Zmass->Fill(mass);
          float pull=fabs(mass-91.1876);
          
          if(pull < minpullZee) {
            minpullZee = pull;
            electrons[0] = iele1;
            electrons[1] = iele2;
          }
        }
      }
    }

    float minpullZmumu = 1000.;
    muons[0]=-1, muons[1]=-1;
    if(nMuon>=2) { 
      for(int imu1=0; imu1<nMuon; imu1++) {
        TLorentzVector muon1(pxMuon[imu1],pyMuon[imu1],pzMuon[imu1],energyMuon[imu1]);
        for(int imu2=imu1+1; imu2<nMuon; imu2++) {
          TLorentzVector muon2(pxMuon[imu2],pyMuon[imu2],pzMuon[imu2],energyMuon[imu2]);
          
          float mass = (muon1+muon2).M();
          float pull=fabs(mass-91.1876);
          
          if(pull < minpullZmumu) {
            minpullZmumu = pull;
            muons[0] = imu1;
            muons[1] = imu2;
          }
        }
      }
    }

    bool Ztag = false;
    // require the reconstruction of the Z
    // the two electrons in the acceptance, pTmin and Z mass
    // or the same with muons
    if (electrons[0]=!-1 && electrons[1]!=-1) {

      TVector3 p3Ele1(pxEle[electrons[0]],pyEle[electrons[0]],pzEle[electrons[0]]);
      TVector3 p3Ele2(pxEle[electrons[1]],pyEle[electrons[1]],pzEle[electrons[1]]);

      if( minpullZee < 9.0 &&
          fabs(etaEle[electrons[0]]) < 2.5 &&
          fabs(etaEle[electrons[1]]) < 2.5 &&
          p3Ele1.Pt() > 5.0 &&
          p3Ele2.Pt() > 5.0 ) Ztag = true;

    } else if (muons[0]!=-1 && muons[1]!=-1) {

      TVector3 p3Mu1(pxEle[muons[0]],pyEle[muons[0]],pzEle[muons[0]]);
      TVector3 p3Mu2(pxEle[muons[1]],pyEle[muons[1]],pzEle[muons[1]]);
      
      if( minpullZmumu < 9.0 &&
          fabs(etaMuon[muons[0]]) < 2.4 &&
          fabs(etaMuon[muons[1]]) < 2.4 &&
          p3Mu1.Pt() > 5.0 &&
          p3Mu2.Pt() > 5.0 ) Ztag = true;

    }


    if(!Ztag) continue;
    
    cout << "Looking for fakes" << endl;

    // find fake electrons
    for(int iele=0;iele<nEle;iele++) {

      cout << "iele = " << iele << "  electrons[0] = " << electrons[0] << "   electrons[1] = " << electrons[1] << endl;

      if( iele==electrons[0] || iele==electrons[1] ) continue;

      TVector3 eleP3(pxEle[iele],pyEle[iele],pzEle[iele]);

      if( fabs(etaEle[iele]) > 2.5 || eleP3.Pt() < 5.0 ) continue;
      
      AllFakesEta->Fill(etaEle[iele]);
      AllFakesPt->Fill(eleP3.Pt());

      // evaluate the purity of the jets: p=#ele(unmatched)/#ele
      if(mcEleMinusIndex>-1 && mcElePlusIndex>-1) {
        float dRPlus = mcElePlus.DeltaR(eleP3);
        float dRMinus = mcEleMinus.DeltaR(eleP3);
        if(dRPlus>0.3 && dRMinus>0.3) {
          UnmatchedFakesEta->Fill(etaEle[iele]);
          UnmatchedFakesPt->Fill(eleP3.Pt());
        }
      }
      
      /// define the bins in which can be splitted the PDFs
      int iecal = (fabs( etaEle[iele])<1.479) ? 0 : 1;
      int iptbin = (eleP3.Pt()<15.0) ? 0 : 1;
      
      int fullclassRaw = classificationEle[iele];
        
      int iclass = -1;
      int ifullclass = -1;
      if ( fullclassRaw == 0 || fullclassRaw == 100 ) { // golden
        iclass = 0;
        ifullclass = 0;
      }
      else if ( fullclassRaw == 10 || fullclassRaw == 110 ) { // bigbrem
        iclass = 0;
        ifullclass = 1;
      }
      else if ( fullclassRaw == 20 || fullclassRaw == 120 ) { // narrow
        iclass = 0;
        ifullclass = 2;
      }
      else if ( (fullclassRaw >= 30 && fullclassRaw <= 40) ||
                (fullclassRaw >= 130 && fullclassRaw <= 140) ) { // showering + cracks
        iclass = 1;
        ifullclass = 3;
      }

      int sc = superClusterIndexEle[iele];
      float sigmaIEtaIEta = sqrt(fabs(covIEtaIEtaSC[sc]));
      float sigmaIEtaIPhi = sqrt(fabs(covIEtaIPhiSC[sc]));
      float sigmaIPhiIPhi = sqrt(fabs(covIPhiIPhiSC[sc]));
      float s1s9 = eMaxSC[sc]/e3x3SC[sc];
      float s9s25 = e3x3SC[sc]/e5x5SC[sc];
      double dPhiCalo = deltaPhiAtCaloEle[iele];
      double dPhiVtx = deltaPhiAtVtxEle[iele];
      double dEta = deltaEtaAtVtxEle[iele];
      double EoPout = eSeedOverPoutEle[iele];
      double EoP = eSuperClusterOverPEle[iele];
      double HoE = hOverEEle[iele];

      dPhiCaloUnsplitEle    [iecal][iptbin] -> Fill ( dPhiCalo );
      dPhiVtxUnsplitEle     [iecal][iptbin] -> Fill ( dPhiVtx );
      dEtaUnsplitEle        [iecal][iptbin] -> Fill ( dEta );
      EoPoutUnsplitEle      [iecal][iptbin] -> Fill ( EoPout );
      HoEUnsplitEle         [iecal][iptbin] -> Fill ( HoE );
      sigmaIEtaIEtaUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiUnsplitEle [iecal][iptbin] -> Fill ( sigmaIPhiIPhi );
      s1s9UnsplitEle        [iecal][iptbin] -> Fill ( s1s9 );
      s9s25UnsplitEle       [iecal][iptbin] -> Fill ( s9s25 );

      dPhiCaloClassEle    [iecal][iptbin][iclass] -> Fill ( dPhiCalo );
      dPhiVtxClassEle     [iecal][iptbin][iclass] -> Fill ( dPhiVtx );
      dEtaClassEle        [iecal][iptbin][iclass] -> Fill ( dEta );
      EoPoutClassEle      [iecal][iptbin][iclass] -> Fill ( EoPout );
      HoEClassEle         [iecal][iptbin][iclass] -> Fill ( HoE );
      sigmaIEtaIEtaClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiClassEle [iecal][iptbin][iclass] -> Fill ( sigmaIPhiIPhi );
      s1s9ClassEle        [iecal][iptbin][iclass] -> Fill ( s1s9 );
      s9s25ClassEle       [iecal][iptbin][iclass] -> Fill ( s9s25 );

      dPhiCaloFullclassEle    [iecal][iptbin][ifullclass] -> Fill ( dPhiCalo );
      dPhiVtxFullclassEle     [iecal][iptbin][ifullclass] -> Fill ( dPhiVtx );
      dEtaFullclassEle        [iecal][iptbin][ifullclass] -> Fill ( dEta );
      EoPoutFullclassEle      [iecal][iptbin][ifullclass] -> Fill ( EoPout );
      HoEFullclassEle         [iecal][iptbin][ifullclass] -> Fill ( HoE );
      sigmaIEtaIEtaFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIEta );
      sigmaIEtaIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIEtaIPhi );
      sigmaIPhiIPhiFullclassEle [iecal][iptbin][ifullclass] -> Fill ( sigmaIPhiIPhi );
      s1s9FullclassEle        [iecal][iptbin][ifullclass] -> Fill ( s1s9 );
      s9s25FullclassEle       [iecal][iptbin][ifullclass] -> Fill ( s9s25 );

    } // loop on electrons

  } // loop on events

  char filename[200];
  sprintf(filename,"%s-JetPurityEta.root",outname);
  EfficiencyEvaluator PurityOfFakesEta(filename);
  PurityOfFakesEta.AddNumerator(AllFakesEta);
  PurityOfFakesEta.AddNumerator(UnmatchedFakesEta);
  PurityOfFakesEta.SetDenominator(AllFakesEta);
  PurityOfFakesEta.ComputeEfficiencies();
  PurityOfFakesEta.Write();

  sprintf(filename,"%s-JetPurityPt.root",outname);
  EfficiencyEvaluator PurityOfFakesPt(filename);
  PurityOfFakesPt.AddNumerator(AllFakesPt);
  PurityOfFakesPt.AddNumerator(UnmatchedFakesPt);
  PurityOfFakesPt.SetDenominator(AllFakesPt);
  PurityOfFakesPt.ComputeEfficiencies();
  PurityOfFakesPt.Write();
  
}


void LHPdfsProducer::bookHistos() {

  m_Zmass = new TH1F("Zmass", "Zmass", 260, 0., 130.);
  
  int nbins = 100;

  float dPhiCaloMin = -0.3;
  float dPhiCaloMax =  0.3;
  float dPhiVtxMin  = -0.1; // pixelMatchGsfElectron pre-selection: |dPhi| < 0.1
  float dPhiVtxMax  =  0.1;
  float dEtaMin     = -0.02;
  float dEtaMax     =  0.02;
  float EoPoutMin   =  0.0;
  float EoPoutMax   =  8.0;
  float HoEMin      =  0.0; // zero-suppression in HCAL
  float HoEMax      =  0.1; // ??
  float sigmaIEtaIEtaMin = 0.0;
  float sigmaIEtaIEtaMax = 0.05;
  float sigmaIEtaIPhiMin = 0.0;
  float sigmaIEtaIPhiMax = 0.05;
  float sigmaIPhiIPhiMin = 0.0;
  float sigmaIPhiIPhiMax = 0.05;
  float s1s9Min  = 0.0;
  float s1s9Max  = 1.0;
  float s9s25Min = 0.5;
  float s9s25Max = 1.0;
//   float dxyMin    = -0.04;
//   float dxyMax    = 0.04;
//   float dxySigMin    = -10.0;
//   float dxySigMax    = 10.0;

  // booking histos eleID
  // iecal = 0 --> barrel
  // iecal = 1 --> endcap
  for (int iecal=0; iecal<2; iecal++) {

    // iptbin = 0: < 15 GeV
    // iptbin = 1: > 15 GeV
    for(int iptbin=0; iptbin<2; iptbin++) {
      
      char histo[200];
      
      sprintf(histo,"dPhiCaloUnsplit_electrons_%d_%d",iecal,iptbin);
      dPhiCaloUnsplitEle[iecal][iptbin]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
      sprintf(histo,"dPhiVtxUnsplit_electrons_%d_%d",iecal,iptbin);
      dPhiVtxUnsplitEle[iecal][iptbin]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
      sprintf(histo,"dEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      dEtaUnsplitEle[iecal][iptbin]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
      sprintf(histo,"EoPoutUnsplit_electrons_%d_%d",iecal,iptbin);
      EoPoutUnsplitEle[iecal][iptbin]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
      sprintf(histo,"HoEUnsplit_electrons_%d_%d",iecal,iptbin);
      HoEUnsplitEle[iecal][iptbin]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
      sprintf(histo,"sigmaIEtaIEtaUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaIEtaIEtaUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaMin, sigmaIEtaIEtaMax);
      sprintf(histo,"sigmaIEtaIPhiUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaIEtaIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIEtaIPhiMin, sigmaIEtaIPhiMax);
      sprintf(histo,"sigmaIPhiIPhiUnsplit_electrons_%d_%d",iecal,iptbin);
      sigmaIPhiIPhiUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiMin, sigmaIPhiIPhiMax);
      sprintf(histo,"s1s9Unsplit_electrons_%d_%d",iecal,iptbin);
      s1s9UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
      sprintf(histo,"s9s25Unsplit_electrons_%d_%d",iecal,iptbin);
      s9s25UnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
//       dxyUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, dxyMin, dxyMax);
//       sprintf(histo,"dxySigUnsplit_electrons_%d_%d",iecal,iptbin);
//       dxySigUnsplitEle[iecal][iptbin] = new TH1F(histo, histo, nbins, dxySigMin, dxySigMax);

      // iclass = 0: non-showering
      // iclass = 1: showering
      for(int iclass=0; iclass<2; iclass++) {
      
	sprintf(histo,"dPhiCaloClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dPhiCaloClassEle[iecal][iptbin][iclass]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
	sprintf(histo,"dPhiVtxClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dPhiVtxClassEle[iecal][iptbin][iclass]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
	sprintf(histo,"dEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	dEtaClassEle[iecal][iptbin][iclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPoutClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	EoPoutClassEle[iecal][iptbin][iclass]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
	sprintf(histo,"HoEClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	HoEClassEle[iecal][iptbin][iclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
	sprintf(histo,"sigmaIEtaIEtaClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaMin, sigmaIEtaIEtaMax);
	sprintf(histo,"sigmaIEtaIPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaIEtaIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIEtaIPhiMin, sigmaIEtaIPhiMax);
	sprintf(histo,"sigmaIPhiIPhiClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
	sigmaIPhiIPhiClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiMin, sigmaIPhiIPhiMax);
	sprintf(histo,"s1s9Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	s1s9ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
	sprintf(histo,"s9s25Class_electrons_%d_%d_%d",iecal,iptbin,iclass);
	s9s25ClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
//         sprintf(histo,"dxyClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
//         dxyClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, dxyMin, dxyMax);
//         sprintf(histo,"dxySigClass_electrons_%d_%d_%d",iecal,iptbin,iclass);
//         dxySigClassEle[iecal][iptbin][iclass] = new TH1F(histo, histo, nbins, dxySigMin, dxySigMax);

      }

      // iclass = 0: golden
      // iclass = 1: bigbrem
      // iclass = 2: narrow
      // iclass = 3: showering + cracks

      for(int ifullclass=0; ifullclass<4; ifullclass++) {
      
	sprintf(histo,"dPhiCaloFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dPhiCaloFullclassEle[iecal][iptbin][ifullclass]    = new TH1F(histo, histo, nbins, dPhiCaloMin, dPhiCaloMax);
	sprintf(histo,"dPhiVtxFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dPhiVtxFullclassEle[iecal][iptbin][ifullclass]     = new TH1F(histo, histo, nbins, dPhiVtxMin, dPhiVtxMax);
	sprintf(histo,"dEtaFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	dEtaFullclassEle[iecal][iptbin][ifullclass]        = new TH1F(histo, histo, nbins, dEtaMin, dEtaMax);
	sprintf(histo,"EoPoutFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	EoPoutFullclassEle[iecal][iptbin][ifullclass]      = new TH1F(histo, histo, nbins, EoPoutMin, EoPoutMax);
	sprintf(histo,"HoEFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	HoEFullclassEle[iecal][iptbin][ifullclass]         = new TH1F(histo, histo, nbins, HoEMin, HoEMax);
	sprintf(histo,"sigmaIEtaIEtaFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaIEtaIEtaFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaIEtaIEtaMin, sigmaIEtaIEtaMax);
	sprintf(histo,"sigmaIEtaIPhiFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaIEtaIPhiFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaIEtaIPhiMin, sigmaIEtaIPhiMax);
	sprintf(histo,"sigmaIPhiIPhiFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	sigmaIPhiIPhiFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, sigmaIPhiIPhiMin, sigmaIPhiIPhiMax);
	sprintf(histo,"s1s9Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	s1s9FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, s1s9Min, s1s9Max);
	sprintf(histo,"s9s25Fullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
	s9s25FullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, s9s25Min, s9s25Max);
//         sprintf(histo,"dxyFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
//         dxyFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, dxyMin, dxyMax);
//         sprintf(histo,"dxySigFullclass_electrons_%d_%d_%d",iecal,iptbin,ifullclass);
//         dxySigFullclassEle[iecal][iptbin][ifullclass] = new TH1F(histo, histo, nbins, dxySigMin, dxySigMax);

      }

    }

  }

}

void LHPdfsProducer::saveHistos(const char *filename) {

  char outfilename[200];
  sprintf(outfilename,"%s",filename);
  TFile *file = TFile::Open(outfilename,"recreate");
  file->mkdir("pdfsProducer","pdfs created from trees");
  file->cd("pdfsProducer");

  m_Zmass->Write();

  for (int iecal=0; iecal<2; iecal++) {

    for(int iptbin=0; iptbin<2; iptbin++) {
      
      dPhiCaloUnsplitEle[iecal][iptbin]->Write();
      dPhiVtxUnsplitEle[iecal][iptbin]->Write();
      dEtaUnsplitEle[iecal][iptbin]->Write();
      EoPoutUnsplitEle[iecal][iptbin]->Write();
      HoEUnsplitEle[iecal][iptbin]->Write();
      sigmaIEtaIEtaUnsplitEle[iecal][iptbin]->Write();
      sigmaIEtaIPhiUnsplitEle[iecal][iptbin]->Write();
      sigmaIPhiIPhiUnsplitEle[iecal][iptbin]->Write();
      s1s9UnsplitEle[iecal][iptbin]->Write();
      s9s25UnsplitEle[iecal][iptbin]->Write();
//       dxyUnsplitEle[iecal][iptbin]->Write();
//       dxySigUnsplitEle[iecal][iptbin]->Write();

      for(int iclass=0; iclass<2; iclass++) {
      
	dPhiCaloClassEle[iecal][iptbin][iclass]->Write();
	dPhiVtxClassEle[iecal][iptbin][iclass]->Write();
	dEtaClassEle[iecal][iptbin][iclass]->Write();
	EoPoutClassEle[iecal][iptbin][iclass]->Write();
	HoEClassEle[iecal][iptbin][iclass]->Write();
	sigmaIEtaIEtaClassEle[iecal][iptbin][iclass]->Write();
	sigmaIEtaIPhiClassEle[iecal][iptbin][iclass]->Write();
	sigmaIPhiIPhiClassEle[iecal][iptbin][iclass]->Write();
	s1s9ClassEle[iecal][iptbin][iclass]->Write();
	s9s25ClassEle[iecal][iptbin][iclass]->Write();
// 	dxyClassEle[iecal][iptbin][iclass]->Write();
// 	dxySigClassEle[iecal][iptbin][iclass]->Write();

      }

      for(int ifullclass=0; ifullclass<4; ifullclass++) {
      
	dPhiCaloFullclassEle[iecal][iptbin][ifullclass]->Write();
	dPhiVtxFullclassEle[iecal][iptbin][ifullclass]->Write();
	dEtaFullclassEle[iecal][iptbin][ifullclass]->Write();
	EoPoutFullclassEle[iecal][iptbin][ifullclass]->Write();
	HoEFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaIEtaIEtaFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaIEtaIPhiFullclassEle[iecal][iptbin][ifullclass]->Write();
	sigmaIPhiIPhiFullclassEle[iecal][iptbin][ifullclass]->Write();
	s1s9FullclassEle[iecal][iptbin][ifullclass]->Write();
	s9s25FullclassEle[iecal][iptbin][ifullclass]->Write();
// 	dxyFullclassEle[iecal][iptbin][ifullclass]->Write();
// 	dxySigFullclassEle[iecal][iptbin][ifullclass]->Write();

      }

    }

  }

  file->Close();
  
}

void LHPdfsProducer::isEleID(int eleIndex, bool *eleIdOutput, bool *isolOutput, bool *convRejOutput) {

  *eleIdOutput = *isolOutput = *convRejOutput = false;

  Utils anaUtils;
  int gsf = gsfTrackIndexEle[eleIndex];
  TVector3 p3Ele(pxEle[eleIndex],pyEle[eleIndex],pzEle[eleIndex]);
  float pt = p3Ele.Pt();

  // if is ECAL driven, take the electron ID variables from the standard electron
  // above all, take the ECAL supercluster instead of PF super cluster
  float HoE, s9s25, deta, dphiin, dphiout, fbrem, see, spp, eopout, eop;
  float e1, e4SwissCross, fidFlagSC, seedRecHitFlag, seedTime, seedChi2;
  bool ecaldriven = anaUtils.electronRecoType(recoFlagsEle[eleIndex], isEcalDriven);
  HoE = hOverEEle[eleIndex];
  deta = deltaEtaAtVtxEle[eleIndex];
  dphiin = deltaPhiAtVtxEle[eleIndex];
  dphiout = deltaPhiAtCaloEle[eleIndex];
  fbrem = fbremEle[eleIndex];
  eopout = eSeedOverPoutEle[eleIndex];
  eop = eSuperClusterOverPEle[eleIndex];
  if(ecaldriven) {
    int sc = superClusterIndexEle[eleIndex];
    s9s25 = e3x3SC[sc]/e5x5SC[sc];
    see = sqrt(covIEtaIEtaSC[sc]);
    spp = sqrt(covIPhiIPhiSC[sc]);
    e1 = eMaxSC[sc];
    e4SwissCross = e4SwissCrossSC[sc];
    fidFlagSC = fiducialFlagsEle[eleIndex];
    seedRecHitFlag = recoFlagSC[sc];
    seedTime = timeSC[sc];
    seedChi2 = chi2SC[sc];
  } else {
    int sc = PFsuperClusterIndexEle[eleIndex];
    if(sc>-1) {
      s9s25 = e3x3PFSC[sc]/e5x5PFSC[sc];
      see = sqrt(covIEtaIEtaPFSC[sc]);
      spp = sqrt(covIPhiIPhiPFSC[sc]);
      e1 = eMaxSC[sc];
      e4SwissCross = e4SwissCrossSC[sc];
      fidFlagSC = fiducialFlagsEle[eleIndex];
      seedRecHitFlag = recoFlagSC[sc];
      seedTime = timeSC[sc];
      seedChi2 = chi2SC[sc];
    } else {
      s9s25 = 999.;
      see = 999.;
      spp = 999.;
    }
  }

  EgammaCutBasedID.SetEcalFiducialRegion( fiducialFlagsEle[eleIndex] );
  EgammaCutBasedID.SetRecoFlag(recoFlagsEle[eleIndex]);
  EgammaCutBasedID.applyElectronIDOnPFlowElectrons(true);
  EgammaCutBasedID.SetHOverE( HoE );
  EgammaCutBasedID.SetS9S25( s9s25 );
  EgammaCutBasedID.SetDEta( deta );
  EgammaCutBasedID.SetDPhiIn( dphiin );
  EgammaCutBasedID.SetDPhiOut( dphiout );
  EgammaCutBasedID.SetBremFraction( fbrem );
  EgammaCutBasedID.SetSigmaEtaEta( see );
  EgammaCutBasedID.SetSigmaPhiPhi( spp );
  EgammaCutBasedID.SetEOverPout( eopout );
  EgammaCutBasedID.SetEOverPin( eop );
  EgammaCutBasedID.SetElectronClass ( classificationEle[eleIndex] );
  //  EgammaCutBasedID.SetEgammaCutBasedID ( anaUtils.electronIdVal(eleIdCutsEle[eleIndex],eleIdLoose) );
  EgammaCutBasedID.SetLikelihood( eleIdLikelihoodEle[eleIndex] );
  EgammaCutBasedID.SetEcalIsolation( dr03EcalRecHitSumEtEle[eleIndex] );
  EgammaCutBasedID.SetTrkIsolation( dr03TkSumPtEle[eleIndex] );
  EgammaCutBasedID.SetHcalIsolation( dr03HcalTowerSumEtEle[eleIndex] );
  EgammaCutBasedID.SetCombinedIsolation( (dr03TkSumPtEle[eleIndex] + 
                                          TMath::Max(0.0,dr03EcalRecHitSumEtEle[eleIndex]-1.0) + 
                                          dr03HcalTowerSumEtEle[eleIndex]) / pt );
  EgammaCutBasedID.SetMissingHits( expInnerLayersGsfTrack[gsf] );
  EgammaCutBasedID.SetConvDist( fabs(convDistEle[eleIndex]) );
  EgammaCutBasedID.SetConvDcot( fabs(convDcotEle[eleIndex]) );

  // ECAL cleaning variables
  EgammaCutBasedID.m_cleaner->SetE1(e1);
  EgammaCutBasedID.m_cleaner->SetE4SwissCross(e4SwissCross);
  EgammaCutBasedID.m_cleaner->SetFiducialFlag(fidFlagSC);
  EgammaCutBasedID.m_cleaner->SetSeedFlag(seedRecHitFlag);
  EgammaCutBasedID.m_cleaner->SetSeedTime(seedTime);
  EgammaCutBasedID.m_cleaner->SetSeedChi2(seedChi2);

  //  return EgammaCutBasedID.output(); // class dependent result
  *eleIdOutput = EgammaCutBasedID.outputNoClassEleId();
  *isolOutput = EgammaCutBasedID.outputNoClassIso();
  *convRejOutput = EgammaCutBasedID.outputNoClassConv();

}

void LHPdfsProducer::configSelection(Selection* selection, Counters* counters) {

  m_selection->addSwitch("isData");
  m_selection->addSwitch("goodRunLS");
  m_selection->addSwitch("requireTriggerSignal");
  m_selection->addSwitch("requireTriggerQCDBack");
  m_selection->addSwitch("applyIsolationOnProbe");
  m_selection->addCut("ptHat");
  m_selection->addCut("meeWindow");
  m_selection->addCut("etaEleAcc");
  m_selection->addCut("ptEleAcc");
  m_selection->addCut("etaJetAcc");
  m_selection->addCut("ptJetAcc");
  m_selection->addCut("jetDeltaPhi");
  m_selection->addCut("jetInvMass");
  m_selection->addCut("met");
  m_selection->addCut("antiIsolTracker");
  m_selection->addCut("antiIsolEcal");
  m_selection->addCut("relSumPtTracks");
  m_selection->addCut("spikeFraction");
  m_selection->addStringParameter("electronIDType");
  m_selection->summary();

  m_counters->SetTitle("EVENT COUNTER");
  m_counters->AddVar("allEvents");
  m_counters->AddVar("pthat");
  m_counters->AddVar("trigger");
  m_counters->AddVar("oneele");
  m_counters->AddVar("onejet");
  m_counters->AddVar("eletot");
  m_counters->AddVar("deltaphi");
  m_counters->AddVar("tagandprobe");
  m_counters->AddVar("met");
  m_counters->AddVar("invmass");
  m_counters->AddVar("trackerNotIsol");
  m_counters->AddVar("ecalNotIsol");
  m_counters->AddVar("fullSelection");
}

void LHPdfsProducer::displayEfficiencies(const char *counterSuffix) {

  std::cout << "----------------------------------------" << std::endl;
  std::cout << "+++ DETAILED EFFICIENCY +++ " << std::endl;

  char namefile[500];
  sprintf(namefile,"%s",counterSuffix);
  
  m_counters->Draw();  
  m_counters->Draw("pthat","allEvents");
  m_counters->Draw("trigger","pthat");
  m_counters->Draw("oneele","trigger");
  m_counters->Draw("onejet","oneele");
  m_counters->Draw("eletot","onejet");
  m_counters->Draw("deltaphi","eletot");
  m_counters->Draw("tagandprobe","deltaphi");
  m_counters->Draw("met","tagandprobe");
  m_counters->Draw("invmass","met");
  m_counters->Draw("trackerNotIsol","invmass");
  m_counters->Draw("ecalNotIsol","trackerNotIsol");
  m_counters->Draw("fullSelection","ecalNotIsol");

  m_counters->Save(namefile,"recreate");
}



void LHPdfsProducer::setRequiredTriggers(const std::vector<std::string>& reqTriggers) {
  requiredTriggers=reqTriggers;
}

bool LHPdfsProducer::hasPassedHLT() {
  Utils anaUtils;
  return anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
}

void LHPdfsProducer::setJsonGoodRunList(const string& jsonFilePath) {

  jsonFile=jsonFilePath;
}

void LHPdfsProducer::fillRunLSMap() {
  
  if (jsonFile == "") {
    std::cout << "Cannot fill RunLSMap. json file not configured" << std::endl;
    return;
  }
  
  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile.c_str());
  if (!jsonFileStream.is_open()) {
    std::cout << "Unable to open file " << jsonFile << std::endl;
    return;
  }
  
  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun) {

    const json::Array& lsSegment = (*itRun).element;
    LSSegments thisRunSegments;
    for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator) {

      json::Array lsSegment=(*lsIterator);
      json::Number lsStart=lsSegment[0];
      json::Number lsEnd=lsSegment[1];
      aLSSegment thisSegment;
      thisSegment.first=lsStart.Value();
      thisSegment.second=lsEnd.Value();
      thisRunSegments.push_back(thisSegment);
    }
    goodRunLS.insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
  }

  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS.size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS.begin(); itR!=goodRunLS.end(); ++itR) {
    std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
    for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
      std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] ";
    std::cout << std::endl;
  }
}

bool LHPdfsProducer::isGoodRunLS() {
  
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);

  if (thisRun == goodRunLS.end()) 
    return false;
  
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg) {
    if ( lumiBlock >= (*iSeg).first && lumiBlock <= (*iSeg).second)
      return true;
  }
  return false;
}

bool LHPdfsProducer::reloadTriggerMask(bool newVersion) {

  if(newVersion) {
    std::vector<int> triggerMask;
    for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
      {
        for(unsigned int i=0; i<nameHLT->size(); i++)
          {
            if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
              {
                triggerMask.push_back( indexHLT[i] ) ;
                break;
              }
          }
      }
    m_requiredTriggers = triggerMask;
  } else {
    TString fileName=((TChain*)fChain)->GetFile()->GetName();
    if ( TString(lastFile) != fileName )
      {
	
	std::cout << "[ReloadTriggerMask]::File has changed reloading trigger mask" << std::endl;
        lastFile = fileName;
	
	TTree *treeCond;
	std::cout << "[ReloadTriggerMask]::Opening " << fileName << std::endl;
        treeCond = (TTree*)((TChain*)fChain)->GetFile()->Get("Conditions");
        int           nHLT_;
	std::vector<std::string>  *nameHLT_;
	std::vector<unsigned int> *indexHLT_;
	
        //To get the pointers for the vectors                                                                                                                                  
        nameHLT_=0;
        indexHLT_=0;
	
        treeCond->SetBranchAddress("nHLT", &nHLT_);
        treeCond->SetBranchAddress("nameHLT", &nameHLT_);
        treeCond->SetBranchAddress("indexHLT", &indexHLT_);
        treeCond->GetEntry(0);
	
	std::vector<int> triggerMask;
        for (std::vector< std::string >::const_iterator fIter=requiredTriggers.begin();fIter!=requiredTriggers.end();++fIter)
          {
            for(unsigned int i=0; i<nameHLT_->size(); i++)
              {
                if( !strcmp ((*fIter).c_str(), nameHLT_->at(i).c_str() ) )
		  {
                    triggerMask.push_back( indexHLT_->at(i) ) ;
                    break;
                  }
              }
          }
        m_requiredTriggers = triggerMask;
        for (int i=0;i<m_requiredTriggers.size();++i)
	  std::cout << "[ReloadTriggerMask]::Requiring bit " << m_requiredTriggers[i] << " " << requiredTriggers[i] << std::endl;
      }
  }
}
