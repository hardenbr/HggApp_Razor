#include "CommonTools/include/Utils.hh"
#include "EgammaAnalysisTools/include/SuperClusterWSelection.hh"

#include "cajun/json/reader.h"
#include "cajun/json/elements.h"

#include <iostream>
#include <fstream>

using namespace std;

SuperClusterWSelection::SuperClusterWSelection(TTree *tree)
  : EgammaBase(tree) {
  // default do not check on mc-truth
  m_signal = all;

  // common kinematic selections
  std::string theConfigDir       = "config/superclustersel/";
  std::string fileCuts     = theConfigDir + "cuts.txt";
  std::string fileSwitches = theConfigDir + "switches.txt";

  _selection = new Selection(fileCuts,fileSwitches);
  _counters = new Counters();
  ConfigSelection(_selection,_counters);

  // To read good run list!
  if (_selection->getSwitch("goodRunLS") && _selection->getSwitch("isData")) {
    std::string goodRunGiasoneFile       = "config/vecbos/json/goodRunLS.json";
    setJsonGoodRunList(goodRunGiasoneFile);
    fillRunLSMap();
  }

  WToENuDecay = 0;
  isData_ = _selection->getSwitch("isData");
  if(isData_) mcevent.SetData(true);

}

SuperClusterWSelection::~SuperClusterWSelection() {
  delete _counters;
  delete _selection;
}

void SuperClusterWSelection::ConfigSelection(Selection* selection, Counters* counters) {
  
  selection->addSwitch("isData");
  selection->addSwitch("mcTruth");
  selection->addSwitch("goodRunLS");
  selection->addSwitch("trigger");
  selection->addCut("eventRange");
  selection->addCut("etaSCAcc");
  selection->addCut("ptSCAcc");
  selection->addCut("HoE");
  selection->addCut("Mt");
  selection->summary();

  counters->SetTitle("EVENT COUNTER");
  counters->AddVar("event");
  counters->AddVar("mcTruth");
  counters->AddVar("trigger");
  counters->AddVar("recoSCs");
  counters->AddVar("accSCs");
  counters->AddVar("idSCs");
  counters->AddVar("Mt");
  counters->AddVar("recoGSF");
  counters->AddVar("fullSelection");

}

void SuperClusterWSelection::Loop() {

  if(fChain == 0) return;

  char namefile[500];
  sprintf(namefile,"%s-SCStudy.root",_prefix);
  fileOut_ = TFile::Open(namefile, "recreate");
  this->createOutTree();

  cout << "Requested signal type = " << m_signal << endl;

  // loop over entries
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Total number of entries in the chain = " << nentries << std::endl;
  int maxEvents = nentries;
  if(_selection->getSwitch("eventRange")) {
    maxEvents = (int) _selection->getUpperCut("eventRange");
    cout << "WARNING! switch eventRange ON! Will run only on the first " << maxEvents << " events!" << endl;
  }
  
  unsigned int lastLumi=0;
  unsigned int lastRun=0;
  
  for (Long64_t jentry=0; jentry<maxEvents;jentry++) {
    
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;
    
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    _counters->IncrVar("event");
    
    //Good Run selection
    if (isData_ && _selection->getSwitch("goodRunLS") && !isGoodRunLS()) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun = runNumber;
        lastLumi = lumiBlock;
        std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }
    
    if (isData_ && _selection->getSwitch("goodRunLS") && ( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }

    int indexMcEleWToENu = -1;
    if ( !isData_ ) { 
      mcevent.LoadDecay(nMc,idMc,mothMc);
      mcevent.LoadMomentum(pMc,energyMc,thetaMc,phiMc);
      // check if the decay is a W->enu prompt
      indexMcEleWToENu = mcevent.indexEleWPrompt();
      WToENuDecay = (indexMcEleWToENu >-1 ) ? 1 : 0;
    }

    // MC truth
    int mctruth = 0;
    if( m_signal == wjets ) mctruth = WToENuDecay;
    else if( m_signal == wother ) mctruth = !WToENuDecay;
    else mctruth = 1;
    if( _selection->getSwitch("mcTruth") && !mctruth ) continue;
    _counters->IncrVar("mcTruth");

    // trigger 
    Utils anaUtils;
    // bool passedHLT = anaUtils.getTriggersOR(m_requiredTriggers, firedTrg);
    bool passedHLT = true;

    if( _selection->getSwitch("trigger") && !passedHLT) continue;
    _counters->IncrVar("trigger");

    if( nSC==0 ) continue;
    _counters->IncrVar("recoSCs");

    // loop over reconstructed ECAL superclusters and choose the accepted ones
    std::vector<int> acceptSCs;
    for(int isc=0; isc<nSC; isc++) {
      TVector3 pClu;
      pClu.SetMagThetaPhi(energySC[isc],thetaSC[isc],phiSC[isc]);
      if( _selection->getSwitch("etaSCAcc") && !_selection->passCut("etaSCAcc",fabs(etaSC[isc])) ) continue;
      if( _selection->getSwitch("ptSCAcc") && !_selection->passCut("ptSCAcc",pClu.Pt()) ) continue;
      acceptSCs.push_back(isc);
    }

    // loop over accepted SCs and choose the hardest one
    int bestSC=-1;
    float maxEt=-1.0;
    for(int i=0; i<(int)acceptSCs.size(); i++) {
      int isc = acceptSCs[i];
      TVector3 pClu;
      pClu.SetMagThetaPhi(energySC[isc],thetaSC[isc],phiSC[isc]);
      if(pClu.Pt()>maxEt) {
        bestSC = isc;
        maxEt=pClu.Pt();
      }
    }

    if( bestSC<0 ) continue;
    _counters->IncrVar("accSCs");
    
    TVector3 pClu;
    pClu.SetMagThetaPhi(energySC[bestSC],thetaSC[bestSC],phiSC[bestSC]);
      
    if( _selection->getSwitch("HoE") && !_selection->passCut("HoE",hOverESC[bestSC]) ) continue;
    _counters->IncrVar("idSCs");

    TVector3 p3PFMet(pxPFMet[0],pyPFMet[0],0.0);
    TVector3 pT3SC(pClu.Px(),pClu.Py(),0.0);
    float WPFmT = sqrt(2 * pClu.Pt() * p3PFMet.Mag() * (1-cos(pT3SC.Angle(p3PFMet))) );   

    if( _selection->getSwitch("Mt") && !_selection->passCut("Mt",WPFmT) ) continue;
    _counters->IncrVar("Mt");
    _counters->IncrVar("fullSelection");
    
    int matchedSC=0;
    // look for the GSF track
    for(int iele=0; iele<nEle; iele++) {
      int eleSc = superClusterIndexEle[iele];
      if( eleSc == bestSC ) {
        matchedSC = 1;
        break;
      }
    }

    if( matchedSC ) _counters->IncrVar("recoGSF");

    // fill the output tree
    myMt = WPFmT;
    myEta = pClu.Eta();
    myPt = pClu.Pt();
    myRecoGSF = matchedSC;

    outTree_->Fill();
    
  }

  outTree_->Write();
  fileOut_->Close();

}

void SuperClusterWSelection::setJsonGoodRunList(const string& jsonFilePath)
{
  jsonFile=jsonFilePath;
}

void SuperClusterWSelection::fillRunLSMap()
{
  
  if (jsonFile == "")
    {
      std::cout << "Cannot fill RunLSMap. json file not configured" << std::endl;
      return;
    }

  std::ifstream jsonFileStream;
  jsonFileStream.open(jsonFile.c_str());
  if (!jsonFileStream.is_open())
    {
      std::cout << "Unable to open file " << jsonFile << std::endl;
      return;
    }

  json::Object elemRootFile;
  json::Reader::Read(elemRootFile, jsonFileStream);

  for (json::Object::const_iterator itRun=elemRootFile.Begin();itRun!=elemRootFile.End();++itRun)
    {
      const json::Array& lsSegment = (*itRun).element;
      LSSegments thisRunSegments; 
      for (json::Array::const_iterator lsIterator=lsSegment.Begin();lsIterator!=lsSegment.End();++lsIterator)
	{
	  json::Array lsSegment=(*lsIterator);
	  json::Number lsStart=lsSegment[0];	   
	  json::Number lsEnd=lsSegment[1];
	  aLSSegment thisSegment;
	  thisSegment.first=lsStart.Value();
	  thisSegment.second=lsEnd.Value();
	  thisRunSegments.push_back(thisSegment);
	  //	   std::pair<int, int> lsSegment=std::pair<int, int>(atoi(,lsIterator[1]); 
	}
      goodRunLS.insert(aRunsLSSegmentsMapElement(atoi((*itRun).name.c_str()),thisRunSegments));
    }


  std::cout << "[GoodRunLSMap]::Good Run LS map filled with " << goodRunLS.size() << " runs" << std::endl;
  for (runsLSSegmentsMap::const_iterator itR=goodRunLS.begin(); itR!=goodRunLS.end(); ++itR)
    {
      std::cout << "[GoodRunLSMap]::Run " << (*itR).first <<  " LS ranges are: ";
      for (LSSegments::const_iterator iSeg=(*itR).second.begin();iSeg!=(*itR).second.end();++iSeg)
	std::cout << "[" << (*iSeg).first << "," << (*iSeg).second << "] "; 
      std::cout << std::endl;
    }
}

bool SuperClusterWSelection::isGoodRunLS()
{
  runsLSSegmentsMap::const_iterator thisRun=goodRunLS.find(runNumber);
  if (thisRun == goodRunLS.end())
    return false;
  //  std::cout << runNumber << " found in the good run map" << std::endl;
  for (LSSegments::const_iterator iSeg=goodRunLS[runNumber].begin();iSeg!=goodRunLS[runNumber].end();++iSeg)
    {
      //      std::cout << "Range is [" << (*iSeg).first << "," << (*iSeg).second << "]" << std::endl;
      if ( lumiBlock >= (*iSeg).first && lumiBlock <= (*iSeg).second)
	return true;
    }
  return false;
}

void SuperClusterWSelection::createOutTree() {

  outTree_ = new TTree("T1","supercluser tree");
  outTree_->Branch("eta", &myEta, "eta/F");
  outTree_->Branch("pt", &myPt, "pt/F");
  outTree_->Branch("mt", &myMt, "mt/F");
  outTree_->Branch("recoGSF", &myRecoGSF, "recoGSF/I");
}

void SuperClusterWSelection::displayEfficiencies() {
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "+++ DETAILED EFFICIENCY +++ " << std::endl;
  
  char namefile[500];
  sprintf(namefile,"%s-SC-Counters.root",_prefix);
  
  _counters->Draw();
  _counters->Draw("mcTruth","event");
  _counters->Draw("trigger","mcTruth");
  _counters->Draw("recoSCs","trigger");
  _counters->Draw("accSCs","recoSCs");
  _counters->Draw("idSCs","accSCs");
  _counters->Draw("Mt","idSCs");
  _counters->Draw("fullSelection","mcTruth");
  _counters->Draw("recoGSF","fullSelection");
  
  _counters->Save(namefile,"recreate");

}
