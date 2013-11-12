//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronIDMVA
//
// Helper Class for applying MVA electron ID selection
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef HIGGSANALYSIS_ElectronIDMVA_H
#define HIGGSANALYSIS_ElectronIDMVA_H

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class ElectronIDMVA {
  public:
    ElectronIDMVA();
    ~ElectronIDMVA(); 

    enum MVAType {
      kBaseline = 0,      // SigmaIEtaIEta, DEtaIn, DPhiIn, FBrem, SigmaIPhiIPhi, NBrem, 
                          // OneOverEMinusOneOverP
      kNoIPInfo,          // kBaseline + EOverP, ESeedClusterOverPout, ESeedClusterOverPIn
      kWithIPInfo,        // kV2 + d0 , IP3d, IP3dSig
      kIDIsoCombined      // new implementation
    };

    void   Initialize(std::string methodName,
                      std::string Subdet0Pt10To20Weights , 
                      std::string Subdet1Pt10To20Weights , 
                      std::string Subdet2Pt10To20Weights,
                      std::string Subdet0Pt20ToInfWeights, 
                      std::string Subdet1Pt20ToInfWeights, 
                      std::string Subdet2Pt20ToInfWeights,
                      ElectronIDMVA::MVAType type );
    Bool_t IsInitialized() const { return fIsInitialized; }
    
    double MVAValue(double ElePt , double EleSCEta,
                    double EleSigmaIEtaIEta,
                    double EleDEtaIn,
                    double EleDPhiIn,
                    double EleHoverE,
                    double EleD0,
                    double EleFBrem,
                    double EleEOverP,
                    double EleESeedClusterOverPout,
                    double EleSigmaIPhiIPhi,
                    double EleNBrem,
                    double EleOneOverEMinusOneOverP,
                    double EleESeedClusterOverPIn,
                    double EleIP3d,
                    double EleIP3dSig );

    double MVAValueWithIso(double ElePt , Double_t EleSCEta,
			   double EleSigmaIEtaIEta,
			   double EleDEtaIn,
			   double EleDPhiIn,
			   double EleHoverE,
			   double EleD0,
			   double EleFBrem,
			   double EleEOverP,
			   double EleESeedClusterOverPout,
			   double EleSigmaIPhiIPhi,
			   double EleNBrem,
			   double EleOneOverEMinusOneOverP,
			   double EleESeedClusterOverPIn,
			   double EleIP3d,
			   double EleIP3dSig,
			   double EleGsfTrackChi2OverNdof,
			   double EledEtaCalo,
			   double EledPhiCalo,
			   double EleR9,
			   double EleSCEtaWidth,
			   double EleSCPhiWidth,
			   double EleCovIEtaIPhi,
			   double ElePreShowerOverRaw,
			   double EleChargedIso03OverPt,
			   double EleNeutralHadronIso03OverPt,
			   double EleGammaIso03OverPt,
			   double EleChargedIso04OverPt,
			   double EleNeutralHadronIso04OverPt,
			   double EleGammaIso04OverPt);
    
    
  protected:
    TMVA::Reader             *fTMVAReader[6];
    std::string               fMethodname;
    MVAType                   fMVAType;
    
    Bool_t                    fIsInitialized;
    
    Float_t                   fMVAVar_EleSigmaIEtaIEta; 
    Float_t                   fMVAVar_EleDEtaIn; 
    Float_t                   fMVAVar_EleDPhiIn; 
    Float_t                   fMVAVar_EleHoverE; 
    Float_t                   fMVAVar_EleD0; 
    Float_t                   fMVAVar_EleFBrem; 
    Float_t                   fMVAVar_EleEOverP; 
    Float_t                   fMVAVar_EleESeedClusterOverPout; 
    Float_t                   fMVAVar_EleSigmaIPhiIPhi; 
    Float_t                   fMVAVar_EleNBrem; 
    Float_t                   fMVAVar_EleOneOverEMinusOneOverP; 
    Float_t                   fMVAVar_EleESeedClusterOverPIn; 
    Float_t                   fMVAVar_EleIP3d; 
    Float_t                   fMVAVar_EleIP3dSig; 
    Float_t                   fMVAVar_EleGsfTrackChi2OverNdof;
    Float_t                   fMVAVar_EledEtaCalo;
    Float_t                   fMVAVar_EledPhiCalo;
    Float_t                   fMVAVar_EleR9;
    Float_t                   fMVAVar_EleSCEtaWidth;
    Float_t                   fMVAVar_EleSCPhiWidth;
    Float_t                   fMVAVar_EleCovIEtaIPhi;
    Float_t                   fMVAVar_ElePreShowerOverRaw;
    Float_t                   fMVAVar_EleChargedIso03OverPt;
    Float_t                   fMVAVar_EleNeutralHadronIso03OverPt;
    Float_t                   fMVAVar_EleGammaIso03OverPt;
    Float_t                   fMVAVar_EleChargedIso04OverPt;
    Float_t                   fMVAVar_EleNeutralHadronIso04OverPt;
    Float_t                   fMVAVar_EleGammaIso04OverPt;
    
};

#endif
