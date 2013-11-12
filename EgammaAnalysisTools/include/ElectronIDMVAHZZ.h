//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronIDMVAHZZ
//
// Helper Class for applying MVA electron ID selection
//
// Authors: D.Benedetti, impl. E.Di Marco
//--------------------------------------------------------------------------------------------------

#ifndef HIGGSANALYSIS_ElectronIDMVAHZZ_H
#define HIGGSANALYSIS_ElectronIDMVAHZZ_H

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class ElectronIDMVAHZZ {
  public:
    ElectronIDMVAHZZ();
    ~ElectronIDMVAHZZ(); 

    enum MVAType {
      kBDTSimpleCat = 0,      // the BDT used in H->ZZ
      kBDTSimpleCatNoIPData, // the BDT used in H->ZZ, trained on DATA 
      kBDTSimpleCatData,     // the BDT used in H->ZZ, including IP, trained on DATA
      kBDTDanV0,             // the new data training, with Daniele's variables 
      kBDTSiV0,              // the new data training, with Si's HWW variables 
      kBDTSiV1,              // the new data training, with Si's HWW 2012 variables
      kBDTSiDanV0,           // the new data training, with Si's HWW 2012 + Daniele's variables
      kBDTSiDanV2,           // the new data training, with Si's HWW 2012 + Daniele's variables, with pruning
      kBDTHWWDanV0,          // the new data training for triggering electrons, with Daniele's variables 
      kBDTHWWSiV0,           // the new data training for triggering electrons, with Si's HWW variables 
      kBDTHWWSiV1,           // the new data training for triggering electrons, with Si's HWW 2012 variables
      kBDTHWWSiDanV0,        // the new data training for triggering electrons, with Si's HWW 2012 + Daniele's variables
      kBDTHWWSiDanV2         // the new data training for triggering electrons, with Si's HWW 2012 + Daniele's variables, with pruning and relaxed preselection for triggering
    };

    void   Initialize(std::string methodName,
                      std::string weights , 
                      ElectronIDMVAHZZ::MVAType type );
    Bool_t IsInitialized() const { return fIsInitialized; }
    
    double MVAValue(double ElePt , double EleSCEta,
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
                    double EleHWWPresel);



  protected:
    TMVA::Reader             *fTMVAReader[1];
    std::string               fMethodname;
    MVAType                   fMVAType;
    
    Bool_t                    fIsInitialized;
    Float_t                   fMVAVar_EleSCEta;
    Float_t                   fMVAVar_ElePt;
    Float_t                   fMVAVar_EleFBrem;     
    Float_t                   fMVAVar_EleNBrems;
    Float_t                   fMVAVar_EleDEtaIn; 
    Float_t                   fMVAVar_EleDPhiIn; 
    Float_t                   fMVAVar_EleDPhiCalo; 
    Float_t                   fMVAVar_EleDEtaCalo; 
    Float_t                   fMVAVar_EleDEtaEleOut;
    Float_t                   fMVAVar_EleSigmaIEtaIEta; 
    Float_t                   fMVAVar_EleSigmaIPhiIPhi; 
    Float_t                   fMVAVar_EleSigmaIEtaIPhi; 
    Float_t                   fMVAVar_EleHoverE; 
    Float_t                   fMVAVar_EleSuperClusterEOverP;
    Float_t                   fMVAVar_EleE1x5E5x5;
    Float_t                   fMVAVar_EleR9;
    Float_t                   fMVAVar_EleESoRaw;
    Float_t                   fMVAVar_EleEseedOverPout;
    Float_t                   fMVAVar_EleEeleOverPout;
    Float_t                   fMVAVar_Ele1oEm1oP;
    Float_t                   fMVAVar_EleKFChi2;
    Float_t                   fMVAVar_EleGSFChi2;
    Float_t                   fMVAVar_EleKFHits;
    Float_t                   fMVAVar_EleMissHits;
    Float_t                   fMVAVar_EleDistConv;
    Float_t                   fMVAVar_EleDcotConv;
    Float_t                   fMVAVar_EleMatchConv;
    Float_t                   fMVAVar_NVtx;
    Float_t                   fMVAVar_EleEcalSeeded;
    Float_t                   fMVAVar_EleEtaWidth;
    Float_t                   fMVAVar_ElePhiWidth;
    Float_t                   fMVAVar_EleD0;
    Float_t                   fMVAVar_EleIP3d;
    Float_t                   fMVAVar_EleIP3dSig;
    Float_t                   fMVAVar_EleDenomFakeSmurf;    
};

#endif
