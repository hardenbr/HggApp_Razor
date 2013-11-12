#include <TFile.h>
#include "EgammaAnalysisTools/include/ElectronIDMVA.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <math.h>
#include <assert.h>

//--------------------------------------------------------------------------------------------------
ElectronIDMVA::ElectronIDMVA() :
fMethodname("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
  for(UInt_t i=0; i<6; ++i) {
    fTMVAReader[i] = 0;
  }
}



//--------------------------------------------------------------------------------------------------
ElectronIDMVA::~ElectronIDMVA()
{
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMVA::Initialize( std::string methodName,
                                std::string Subdet0Pt10To20Weights , 
                                std::string Subdet1Pt10To20Weights , 
                                std::string Subdet2Pt10To20Weights,
                                std::string Subdet0Pt20ToInfWeights,
                                std::string Subdet1Pt20ToInfWeights, 
                                std::string Subdet2Pt20ToInfWeights,
                                ElectronIDMVA::MVAType type) {

  fIsInitialized = kTRUE;
  fMVAType = type;

  fMethodname = methodName;
    
  for(UInt_t i=0; i<6; ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];

    fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
    fTMVAReader[i]->SetVerbose(kTRUE);

    if (type == kBaseline) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                 );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi         );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                 );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP );      
    }
    
    if (type == kNoIPInfo) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                 );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout  );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi         );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                 );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn   );
    }
    if (type == kWithIPInfo) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
      fTMVAReader[i]->AddVariable( "D0",                    &fMVAVar_EleD0                    );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                 );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout  );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi         );
      fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                 );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn   );
      fTMVAReader[i]->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                  );
      fTMVAReader[i]->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig               );
    }

    if (type == kIDIsoCombined) {
      fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      fTMVAReader[i]->AddVariable( "D0",                    &fMVAVar_EleD0                       );
      fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
      fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
      fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
      fTMVAReader[i]->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                     );
      fTMVAReader[i]->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig                  );

      fTMVAReader[i]->AddVariable( "GsfTrackChi2OverNdof",  &fMVAVar_EleGsfTrackChi2OverNdof     );
      fTMVAReader[i]->AddVariable( "dEtaCalo",              &fMVAVar_EledEtaCalo                 );
      fTMVAReader[i]->AddVariable( "dPhiCalo",              &fMVAVar_EledPhiCalo                 );
      fTMVAReader[i]->AddVariable( "R9",                    &fMVAVar_EleR9                       );
      fTMVAReader[i]->AddVariable( "SCEtaWidth",            &fMVAVar_EleSCEtaWidth               );
      fTMVAReader[i]->AddVariable( "SCPhiWidth",            &fMVAVar_EleSCPhiWidth               );
      fTMVAReader[i]->AddVariable( "CovIEtaIPhi",           &fMVAVar_EleCovIEtaIPhi              );
      if (i == 2 || i == 5) {
        fTMVAReader[i]->AddVariable( "PreShowerOverRaw",      &fMVAVar_ElePreShowerOverRaw       );
      }
      fTMVAReader[i]->AddVariable( "ChargedIso03",          &fMVAVar_EleChargedIso03OverPt       );
      fTMVAReader[i]->AddVariable( "NeutralHadronIso03",    &fMVAVar_EleNeutralHadronIso03OverPt );
      fTMVAReader[i]->AddVariable( "GammaIso03",            &fMVAVar_EleGammaIso03OverPt         );
      fTMVAReader[i]->AddVariable( "ChargedIso04",          &fMVAVar_EleChargedIso04OverPt       );
      fTMVAReader[i]->AddVariable( "NeutralHadronIso04",    &fMVAVar_EleNeutralHadronIso04OverPt );
      fTMVAReader[i]->AddVariable( "GammaIso04",            &fMVAVar_EleGammaIso04OverPt         );
    }

    if (i==0) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt10To20Weights );
    if (i==1) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt10To20Weights );
    if (i==2) fTMVAReader[i]->BookMVA(fMethodname , Subdet2Pt10To20Weights );
    if (i==3) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt20ToInfWeights );
    if (i==4) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt20ToInfWeights );
    if (i==5) fTMVAReader[i]->BookMVA(fMethodname , Subdet2Pt20ToInfWeights );

  }

  std::cout << "Electron ID MVA Initialization\n";
  std::cout << "MethodName : " << fMethodname << " , type == " << type << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet2Pt10To20Weights << std::endl;
  std::cout << "Load weights file : " << Subdet0Pt20ToInfWeights << std::endl;
  std::cout << "Load weights file : " << Subdet1Pt20ToInfWeights << std::endl;
  std::cout << "Load weights file : " << Subdet2Pt20ToInfWeights << std::endl;

}


//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(Double_t ElePt , Double_t EleSCEta,
                                 Double_t EleSigmaIEtaIEta,
                                 Double_t EleDEtaIn,
                                 Double_t EleDPhiIn,
                                 Double_t EleHoverE,
                                 Double_t EleD0,
                                 Double_t EleFBrem,
                                 Double_t EleEOverP,
                                 Double_t EleESeedClusterOverPout,
                                 Double_t EleSigmaIPhiIPhi,
                                 Double_t EleNBrem,
                                 Double_t EleOneOverEMinusOneOverP,
                                 Double_t EleESeedClusterOverPIn,
                                 Double_t EleIP3d,
                                 Double_t EleIP3dSig
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (fabs(EleSCEta) < 1.0) subdet = 0;
  else if (fabs(EleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ElePt > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
  fMVAVar_EleDEtaIn = EleDEtaIn;
  fMVAVar_EleDPhiIn = EleDPhiIn;
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleFBrem = EleFBrem;
  fMVAVar_EleEOverP = EleEOverP;
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
  fMVAVar_EleNBrem = EleNBrem;
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
  fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

//   std::cout << "********************************\n";
//   std::cout << "Electron MVA\n";
//   std::cout << ElePt << " no eleeta " << " non ho phi " <<
//     " : " << EleSCEta << " : " << MVABin << std::endl;
//   std::cout << fMVAVar_EleSigmaIEtaIEta << "\n"
//             << fMVAVar_EleDEtaIn << "\n"
//             << fMVAVar_EleDPhiIn << "\n"
//             << fMVAVar_EleHoverE << "\n"
//             << fMVAVar_EleD0 << "\n"
//             << fMVAVar_EleFBrem << "\n"
//             << fMVAVar_EleEOverP << "\n"
//             << fMVAVar_EleESeedClusterOverPout << "\n"
//             << fMVAVar_EleSigmaIPhiIPhi << "\n"
//             << fMVAVar_EleNBrem << "\n"
//             << fMVAVar_EleOneOverEMinusOneOverP << "\n"
//             << fMVAVar_EleESeedClusterOverPIn << "\n"
//             << fMVAVar_EleIP3d << "\n"
//             << fMVAVar_EleIP3dSig << "\n"
//             << std::endl;
//   std::cout << "MVA: " << mva << std::endl;
//   std::cout << "********************************\n";

  return mva;
}


//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValueWithIso(Double_t ElePt , Double_t EleSCEta,
					Double_t EleSigmaIEtaIEta,
					Double_t EleDEtaIn,
					Double_t EleDPhiIn,
					Double_t EleHoverE,
					Double_t EleD0,
					Double_t EleFBrem,
					Double_t EleEOverP,
					Double_t EleESeedClusterOverPout,
					Double_t EleSigmaIPhiIPhi,
					Double_t EleNBrem,
					Double_t EleOneOverEMinusOneOverP,
					Double_t EleESeedClusterOverPIn,
					Double_t EleIP3d,
					Double_t EleIP3dSig,
					Double_t EleGsfTrackChi2OverNdof,
					Double_t EledEtaCalo,
					Double_t EledPhiCalo,
					Double_t EleR9,
					Double_t EleSCEtaWidth,
					Double_t EleSCPhiWidth,
					Double_t EleCovIEtaIPhi,
					Double_t ElePreShowerOverRaw,
					Double_t EleChargedIso03OverPt,
					Double_t EleNeutralHadronIso03OverPt,
					Double_t EleGammaIso03OverPt,
					Double_t EleChargedIso04OverPt,
					Double_t EleNeutralHadronIso04OverPt,
					Double_t EleGammaIso04OverPt) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Int_t subdet = 0;
  if (fabs(EleSCEta) < 1.0) subdet = 0;
  else if (fabs(EleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ElePt > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
  fMVAVar_EleDEtaIn = EleDEtaIn;
  fMVAVar_EleDPhiIn = EleDPhiIn;
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleFBrem = EleFBrem;
  fMVAVar_EleEOverP = EleEOverP;
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
  fMVAVar_EleNBrem = EleNBrem;
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
  fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;
  fMVAVar_EleGsfTrackChi2OverNdof = EleGsfTrackChi2OverNdof;
  fMVAVar_EledEtaCalo = EledEtaCalo;
  fMVAVar_EledPhiCalo = EledPhiCalo;
  fMVAVar_EleR9 = EleR9;
  fMVAVar_EleSCEtaWidth = EleSCEtaWidth;
  fMVAVar_EleSCPhiWidth = EleSCPhiWidth;
  fMVAVar_EleCovIEtaIPhi = EleCovIEtaIPhi;
  fMVAVar_ElePreShowerOverRaw = ElePreShowerOverRaw;
  fMVAVar_EleChargedIso03OverPt = EleChargedIso03OverPt;
  fMVAVar_EleNeutralHadronIso03OverPt = EleNeutralHadronIso03OverPt;
  fMVAVar_EleGammaIso03OverPt = EleGammaIso03OverPt;
  fMVAVar_EleChargedIso04OverPt = EleChargedIso04OverPt;
  fMVAVar_EleNeutralHadronIso04OverPt = EleNeutralHadronIso04OverPt;
  fMVAVar_EleGammaIso04OverPt = EleGammaIso04OverPt;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;
  assert(MVABin >= 0 && MVABin <= 5);
  reader = fTMVAReader[MVABin];
                                                
  mva = reader->EvaluateMVA( fMethodname );

//   std::cout << "********************************\n";
//   std::cout << "Electron MVA\n";
//   std::cout << ElePt << " no eleeta " << " non ho phi " <<
//     " : " << EleSCEta << " : " << MVABin << std::endl;
//   std::cout << fMVAVar_EleSigmaIEtaIEta << "\n"
//             << fMVAVar_EleDEtaIn << "\n"
//             << fMVAVar_EleDPhiIn << "\n"
//             << fMVAVar_EleHoverE << "\n"
//             << fMVAVar_EleD0 << "\n"
//             << fMVAVar_EleFBrem << "\n"
//             << fMVAVar_EleEOverP << "\n"
//             << fMVAVar_EleESeedClusterOverPout << "\n"
//             << fMVAVar_EleSigmaIPhiIPhi << "\n"
//             << fMVAVar_EleNBrem << "\n"
//             << fMVAVar_EleOneOverEMinusOneOverP << "\n"
//             << fMVAVar_EleESeedClusterOverPIn << "\n"
//             << fMVAVar_EleIP3d << "\n"
//             << fMVAVar_EleIP3dSig << "\n"
//             << std::endl;
//   std::cout << "MVA: " << mva << std::endl;
//   std::cout << "********************************\n";

  return mva;
}
