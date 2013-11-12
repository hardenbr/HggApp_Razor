#include <TFile.h>
#include "EgammaAnalysisTools/include/ElectronIDMVAHZZ.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <math.h>
#include <assert.h>

//--------------------------------------------------------------------------------------------------
ElectronIDMVAHZZ::ElectronIDMVAHZZ() :
fMethodname("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
  for(UInt_t i=0; i<1; ++i) {
    fTMVAReader[i] = 0;
  }
}



//--------------------------------------------------------------------------------------------------
ElectronIDMVAHZZ::~ElectronIDMVAHZZ()
{
  for(UInt_t i=0; i<1; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMVAHZZ::Initialize( std::string methodName,
                                   std::string weights , 
                                   ElectronIDMVAHZZ::MVAType type) {

  fIsInitialized = kTRUE;
  fMVAType = type;

  fMethodname = methodName;
    
  for(UInt_t i=0; i<1; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];

    fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
    fTMVAReader[i]->SetVerbose(kTRUE);

    if (type == kBDTSimpleCat) {
      fTMVAReader[i]->AddVariable("fbrem",       &fMVAVar_EleFBrem);
      fTMVAReader[i]->AddVariable("detain",      &fMVAVar_EleDEtaIn);
      fTMVAReader[i]->AddVariable("dphiin",      &fMVAVar_EleDPhiIn);
      fTMVAReader[i]->AddVariable("sieie",       &fMVAVar_EleSigmaIEtaIEta);
      fTMVAReader[i]->AddVariable("hoe",         &fMVAVar_EleHoverE);
      fTMVAReader[i]->AddVariable("eop",         &fMVAVar_EleSuperClusterEOverP);
      fTMVAReader[i]->AddVariable("e1x5e5x5",    &fMVAVar_EleE1x5E5x5);
      fTMVAReader[i]->AddVariable("eleopout",    &fMVAVar_EleEeleOverPout);
      fTMVAReader[i]->AddVariable("detaeleout",  &fMVAVar_EleDEtaEleOut);
      fTMVAReader[i]->AddVariable("kfchi2",      &fMVAVar_EleKFChi2);
      fTMVAReader[i]->AddVariable("kfhits",      &fMVAVar_EleKFHits);
      fTMVAReader[i]->AddVariable("mishits",     &fMVAVar_EleMissHits);
      fTMVAReader[i]->AddVariable("dist",        &fMVAVar_EleDistConv);
      fTMVAReader[i]->AddVariable("dcot",        &fMVAVar_EleDcotConv);
      fTMVAReader[i]->AddVariable("nvtx",        &fMVAVar_NVtx);  // for the new weight file
      fTMVAReader[i]->AddSpectator("eta",        &fMVAVar_EleSCEta);
      fTMVAReader[i]->AddSpectator("pt",         &fMVAVar_ElePt);
      fTMVAReader[i]->AddSpectator("ecalseed",   &fMVAVar_EleEcalSeeded);
    }

    if (type == kBDTSimpleCatNoIPData) {
      fTMVAReader[i]->AddVariable("fbrem",       &fMVAVar_EleFBrem);
      fTMVAReader[i]->AddVariable("detain",      &fMVAVar_EleDEtaIn);
      fTMVAReader[i]->AddVariable("dphiin",      &fMVAVar_EleDPhiIn);
      fTMVAReader[i]->AddVariable("sieie",       &fMVAVar_EleSigmaIEtaIEta);
      fTMVAReader[i]->AddVariable("sc_phiwidth", &fMVAVar_ElePhiWidth);
      fTMVAReader[i]->AddVariable("hoe",         &fMVAVar_EleHoverE);
      fTMVAReader[i]->AddVariable("eop",         &fMVAVar_EleSuperClusterEOverP);
      fTMVAReader[i]->AddVariable("e1x5e5x5",    &fMVAVar_EleE1x5E5x5);
      fTMVAReader[i]->AddVariable("eleopout",    &fMVAVar_EleEeleOverPout);
      fTMVAReader[i]->AddVariable("detaeleout",  &fMVAVar_EleDEtaEleOut);
      fTMVAReader[i]->AddVariable("kfchi2",      &fMVAVar_EleKFChi2);
      fTMVAReader[i]->AddVariable("kfhits",      &fMVAVar_EleKFHits);
      fTMVAReader[i]->AddVariable("dist",        &fMVAVar_EleDistConv);
      fTMVAReader[i]->AddVariable("dcot",        &fMVAVar_EleDcotConv);
      fTMVAReader[i]->AddVariable("nvtx",        &fMVAVar_NVtx);  // for the new weight file
      fTMVAReader[i]->AddSpectator("eta",        &fMVAVar_EleSCEta);
      fTMVAReader[i]->AddSpectator("pt",         &fMVAVar_ElePt);
      fTMVAReader[i]->AddSpectator("ecalseed",   &fMVAVar_EleEcalSeeded);
    }

    if (type == kBDTSimpleCatData) {
      fTMVAReader[i]->AddVariable("fbrem",       &fMVAVar_EleFBrem);
      fTMVAReader[i]->AddVariable("detain",      &fMVAVar_EleDEtaIn);
      fTMVAReader[i]->AddVariable("dphiin",      &fMVAVar_EleDPhiIn);
      fTMVAReader[i]->AddVariable("sieie",       &fMVAVar_EleSigmaIEtaIEta);
      fTMVAReader[i]->AddVariable("sc_etawidth", &fMVAVar_EleEtaWidth);
      fTMVAReader[i]->AddVariable("sc_phiwidth", &fMVAVar_ElePhiWidth);
      fTMVAReader[i]->AddVariable("hoe",         &fMVAVar_EleHoverE);
      fTMVAReader[i]->AddVariable("eop",         &fMVAVar_EleSuperClusterEOverP);
      fTMVAReader[i]->AddVariable("e1x5e5x5",    &fMVAVar_EleE1x5E5x5);
      fTMVAReader[i]->AddVariable("eleopout",    &fMVAVar_EleEeleOverPout);
      fTMVAReader[i]->AddVariable("detaeleout",  &fMVAVar_EleDEtaEleOut);
      fTMVAReader[i]->AddVariable("kfchi2",      &fMVAVar_EleKFChi2);
      fTMVAReader[i]->AddVariable("kfhits",      &fMVAVar_EleKFHits);
      fTMVAReader[i]->AddVariable("dist",        &fMVAVar_EleDistConv);
      fTMVAReader[i]->AddVariable("dcot",        &fMVAVar_EleDcotConv);
      fTMVAReader[i]->AddVariable("d0",          &fMVAVar_EleD0); 
      fTMVAReader[i]->AddVariable("ip3d",        &fMVAVar_EleIP3d); 
      fTMVAReader[i]->AddVariable("ip3ds",       &fMVAVar_EleIP3dSig); 
      fTMVAReader[i]->AddVariable("nvtx",        &fMVAVar_NVtx);  // for the new weight file
      fTMVAReader[i]->AddSpectator("eta",        &fMVAVar_EleSCEta);
      fTMVAReader[i]->AddSpectator("pt",         &fMVAVar_ElePt);
      fTMVAReader[i]->AddSpectator("ecalseed",   &fMVAVar_EleEcalSeeded);
    }

    if(type >= kBDTDanV0 && type <= kBDTSiDanV2) {
      fTMVAReader[i]->AddVariable("fbrem",            &fMVAVar_EleFBrem);
      fTMVAReader[i]->AddVariable("deta",             &fMVAVar_EleDEtaIn);
      fTMVAReader[i]->AddVariable("dphi",             &fMVAVar_EleDPhiIn);
      fTMVAReader[i]->AddVariable("see",              &fMVAVar_EleSigmaIEtaIEta);
      fTMVAReader[i]->AddVariable("etawidth",         &fMVAVar_EleEtaWidth);
      fTMVAReader[i]->AddVariable("phiwidth",         &fMVAVar_ElePhiWidth);
      fTMVAReader[i]->AddVariable("HoE",              &fMVAVar_EleHoverE);
      fTMVAReader[i]->AddVariable("EoP",              &fMVAVar_EleSuperClusterEOverP);
      fTMVAReader[i]->AddVariable("e1x5e5x5",         &fMVAVar_EleE1x5E5x5);
      fTMVAReader[i]->AddVariable("EoPout",           &fMVAVar_EleEseedOverPout);
      fTMVAReader[i]->AddVariable("eleEoPout",        &fMVAVar_EleEeleOverPout); 
      fTMVAReader[i]->AddVariable("detacalo",         &fMVAVar_EleDEtaCalo);
      fTMVAReader[i]->AddVariable("kfchi2",           &fMVAVar_EleKFChi2);
      fTMVAReader[i]->AddVariable("kfhits",           &fMVAVar_EleKFHits);
      fTMVAReader[i]->AddVariable("spp",              &fMVAVar_EleSigmaIPhiIPhi);
      fTMVAReader[i]->AddVariable("IoEmIoP",          &fMVAVar_Ele1oEm1oP);
      fTMVAReader[i]->AddVariable("nbrems",           &fMVAVar_EleNBrems);
      fTMVAReader[i]->AddVariable("R9",               &fMVAVar_EleR9);
      fTMVAReader[i]->AddVariable("dphicalo",         &fMVAVar_EleDPhiCalo);
      fTMVAReader[i]->AddVariable("gsfchi2",          &fMVAVar_EleGSFChi2);
      fTMVAReader[i]->AddVariable("PreShowerOverRaw", &fMVAVar_EleESoRaw);

      fTMVAReader[i]->AddSpectator("eta",             &fMVAVar_EleSCEta);
      fTMVAReader[i]->AddSpectator("pt",              &fMVAVar_ElePt);
      fTMVAReader[i]->AddSpectator("matchConv",       &fMVAVar_EleMatchConv);
    }

    if(type >= kBDTHWWDanV0 && type <= kBDTHWWSiDanV2) {
      fTMVAReader[i]->AddVariable("fbrem",            &fMVAVar_EleFBrem);
      fTMVAReader[i]->AddVariable("deta",             &fMVAVar_EleDEtaIn);
      fTMVAReader[i]->AddVariable("dphi",             &fMVAVar_EleDPhiIn);
      fTMVAReader[i]->AddVariable("see",              &fMVAVar_EleSigmaIEtaIEta);
      fTMVAReader[i]->AddVariable("etawidth",         &fMVAVar_EleEtaWidth);
      fTMVAReader[i]->AddVariable("phiwidth",         &fMVAVar_ElePhiWidth);
      fTMVAReader[i]->AddVariable("HoE",              &fMVAVar_EleHoverE);
      fTMVAReader[i]->AddVariable("EoP",              &fMVAVar_EleSuperClusterEOverP);
      fTMVAReader[i]->AddVariable("e1x5e5x5",         &fMVAVar_EleE1x5E5x5);
      fTMVAReader[i]->AddVariable("EoPout",           &fMVAVar_EleEseedOverPout);
      fTMVAReader[i]->AddVariable("eleEoPout",        &fMVAVar_EleEeleOverPout); 
      fTMVAReader[i]->AddVariable("detacalo",         &fMVAVar_EleDEtaCalo);
      fTMVAReader[i]->AddVariable("kfchi2",           &fMVAVar_EleKFChi2);
      fTMVAReader[i]->AddVariable("kfhits",           &fMVAVar_EleKFHits);
      fTMVAReader[i]->AddVariable("spp",              &fMVAVar_EleSigmaIPhiIPhi);
      fTMVAReader[i]->AddVariable("IoEmIoP",          &fMVAVar_Ele1oEm1oP);
      fTMVAReader[i]->AddVariable("nbrems",           &fMVAVar_EleNBrems);
      fTMVAReader[i]->AddVariable("R9",               &fMVAVar_EleR9);
      fTMVAReader[i]->AddVariable("dphicalo",         &fMVAVar_EleDPhiCalo);
      fTMVAReader[i]->AddVariable("gsfchi2",          &fMVAVar_EleGSFChi2);
      fTMVAReader[i]->AddVariable("PreShowerOverRaw", &fMVAVar_EleESoRaw);
      if(type == kBDTHWWSiDanV2) {
        fTMVAReader[i]->AddVariable("d0",             &fMVAVar_EleD0); 
        fTMVAReader[i]->AddVariable("ip3d",           &fMVAVar_EleIP3d); 
      }

      fTMVAReader[i]->AddSpectator("eta",             &fMVAVar_EleSCEta);
      fTMVAReader[i]->AddSpectator("pt",              &fMVAVar_ElePt);
      if(type != kBDTHWWSiDanV2) fTMVAReader[i]->AddSpectator("DenomFakeSmurf",  &fMVAVar_EleDenomFakeSmurf);
    }

    if (i==0) fTMVAReader[i]->BookMVA(fMethodname , weights );

  }

  std::cout << "Electron ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodname << " , type == " << type << std::endl;
  std::cout << "Load weights file : " << weights << std::endl;

}


//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVAHZZ::MVAValue(double ElePt , double EleSCEta,
                                    double EleFBrem,
                                    double EleNBrems,
                                    double EleDEtaIn,
                                    double EleDPhiIn,
                                    double EleDPhiCalo,
                                    double EleDEtaCalo,
                                    double EleDEtaEleOut,
                                    double EleSigmaIEtaIEta,
                                    double EleSigmaIPhiIPhi,
                                    double EleSigmaIEtaIPhi,
                                    double EleHoverE,
                                    double EleSuperClusterEOverP,
                                    double EleE1x5E5x5,
                                    double EleR9,
                                    double EleESoRaw,
                                    double EleEseedOverPout,
                                    double EleEeleOverPout,
                                    double Ele1oEm1oP,
                                    double EleKFChi2,
                                    double EleKFHits,
                                    double EleGSFChi2,
                                    double EleMissHits,
                                    double EleDistConv,
                                    double EleDcotConv,
                                    double EleMatchConv,
                                    double NVtx,
                                    double EleEcalSeeded,
                                    double EleEtaWidth,
                                    double ElePhiWidth,
                                    double EleD0,
                                    double EleIP3d,
                                    double EleIP3dSig,
                                    double EleHWWPresel) {


  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVAHZZ not properly initialized.\n"; 
    return -9999;
  }

  //set all input variables
  fMVAVar_EleSCEta = EleSCEta;
  fMVAVar_ElePt = ElePt;
  fMVAVar_EleFBrem = EleFBrem;     
  fMVAVar_EleNBrems = EleNBrems;
  fMVAVar_EleDEtaIn = EleDEtaIn; 
  fMVAVar_EleDPhiIn = EleDPhiIn; 
  fMVAVar_EleDPhiCalo = EleDPhiCalo; 
  fMVAVar_EleDEtaCalo = EleDEtaCalo; 
  fMVAVar_EleDEtaEleOut = EleDEtaEleOut;
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta; 
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi; 
  fMVAVar_EleSigmaIEtaIPhi = EleSigmaIEtaIPhi; 
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleSuperClusterEOverP = EleSuperClusterEOverP;
  fMVAVar_EleE1x5E5x5 = EleE1x5E5x5;
  fMVAVar_EleR9 = EleR9;
  fMVAVar_EleESoRaw = EleESoRaw;
  fMVAVar_EleEseedOverPout = EleEseedOverPout;
  fMVAVar_EleEeleOverPout = EleEeleOverPout;
  fMVAVar_Ele1oEm1oP = Ele1oEm1oP;
  fMVAVar_EleKFChi2 = EleKFChi2;
  fMVAVar_EleGSFChi2 = EleGSFChi2;
  fMVAVar_EleKFHits = EleKFHits;
  fMVAVar_EleMissHits = EleMissHits;
  fMVAVar_EleDistConv = EleDistConv;
  fMVAVar_EleDcotConv = EleDcotConv;
  fMVAVar_EleMatchConv = EleMatchConv;
  fMVAVar_NVtx = NVtx;
  fMVAVar_EleEcalSeeded = EleEcalSeeded;
  fMVAVar_EleEtaWidth = EleEtaWidth;
  fMVAVar_ElePhiWidth = ElePhiWidth;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;
  fMVAVar_EleDenomFakeSmurf = EleHWWPresel;

  // apply the boundaries used in the training
  if(fMVAVar_EleFBrem < -1.)
    fMVAVar_EleFBrem = -1.;

  fMVAVar_EleDEtaIn = fabs(fMVAVar_EleDEtaIn);
  if(fMVAVar_EleDEtaIn > 0.06)
    fMVAVar_EleDEtaIn = 0.06;

  fMVAVar_EleDPhiIn = fabs(fMVAVar_EleDPhiIn);
  if(fMVAVar_EleDPhiIn > 0.6)
    fMVAVar_EleDPhiIn = 0.6;

  fMVAVar_EleDEtaCalo = fabs(fMVAVar_EleDEtaCalo);
  if(fMVAVar_EleDEtaCalo > 0.2)
    fMVAVar_EleDEtaCalo = 0.2;

  fMVAVar_EleDPhiCalo = fabs(fMVAVar_EleDPhiCalo);
  if(fMVAVar_EleDPhiCalo > 0.4)
    fMVAVar_EleDPhiCalo = 0.4;

  if(fMVAVar_EleEeleOverPout > 20.)
    fMVAVar_EleEeleOverPout = 20;

  if(fMVAVar_EleEseedOverPout > 20.)
    fMVAVar_EleEseedOverPout = 20;

  if(fMVAVar_EleSuperClusterEOverP > 20.)
    fMVAVar_EleSuperClusterEOverP = 20.;

  fMVAVar_EleDEtaEleOut = fabs(fMVAVar_EleDEtaEleOut);
  if(fMVAVar_EleDEtaEleOut > 0.2)
    fMVAVar_EleDEtaEleOut = 0.2;

  if(fMVAVar_EleKFChi2 < 0.)
    fMVAVar_EleKFChi2 = 0.;

  if(fMVAVar_EleKFChi2 > 10.)
    fMVAVar_EleKFChi2 = 10.;

  if(fMVAVar_EleGSFChi2 > 200.)
    fMVAVar_EleGSFChi2 = 200.;

  if(fMVAVar_EleE1x5E5x5 < -1.)
    fMVAVar_EleE1x5E5x5 = -1;

  if(fMVAVar_EleE1x5E5x5 > 2.)
    fMVAVar_EleE1x5E5x5 = 2.;

  if(fMVAVar_EleR9 > 5.)
    fMVAVar_EleR9 = 5.;

  if(fMVAVar_EleDistConv > 15.)
    fMVAVar_EleDistConv = 15.;
  if(fMVAVar_EleDistConv < -15.)
    fMVAVar_EleDistConv = -15.;
  fMVAVar_EleDistConv = fabs(fMVAVar_EleDistConv);

  if(fMVAVar_EleDcotConv > 3.)
    fMVAVar_EleDcotConv = 3.;
  if(fMVAVar_EleDcotConv < -3.)
    fMVAVar_EleDcotConv = -3.;
  fMVAVar_EleDcotConv = fabs(fMVAVar_EleDcotConv);

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = 0;
  // assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
  
                                                
  mva = reader->EvaluateMVA( fMethodname );

//   std::cout << "********************************\n";
//   std::cout << "Electron MVA\n";
//   std::cout << ElePt << " no eleeta " << " non ho phi " <<
//     " : " << EleSCEta << " : " << MVABin << std::endl;
//   std::cout << fMVAVar_EleFBrem << endl     
//   << fMVAVar_EleDEtaIn << endl 
//   << fMVAVar_EleDPhiIn << endl 
//   << fMVAVar_EleDEtaEleOut << endl
//   << fMVAVar_EleSigmaIEtaIEta << endl 
//   << fMVAVar_EleHoverE << endl 
//   << fMVAVar_EleE1x5E5x5 << endl
//   << fMVAVar_EleEOverPout << endl
//   << fMVAVar_EleKFChi2 << endl
//   << fMVAVar_EleKFHits << endl
//   << fMVAVar_EleMissHits << endl
//   << fMVAVar_EleDistConv << endl
//   << fMVAVar_EleDcotConv << endl
//   << fMVAVar_NVtx << endl
//   << fMVAVar_EleEcalSeeded << endl
//   << std::endl;
//   std::cout << "MVA: " << mva << std::endl;
//   std::cout << "********************************\n";

  return mva;
}
