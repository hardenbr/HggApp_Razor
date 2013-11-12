#ifndef ELECTRONBESTCANDIDATESELECTOR_HH
#define ELECTRONBESTCANDIDATESELECTOR_HH

#include <vector>
#include <functional>

struct ElectronQualityData {
  
  int index;
  float SCenergy;
  float ptAtVtx;
  float trackerSumPt;
  float ecalSumEt;
  float hcalSumEt;
  float electronIdLH;

};

class ElectronBestCandidateSelector {
  
public:

  ElectronBestCandidateSelector( std::vector<ElectronQualityData> );
  virtual ~ElectronBestCandidateSelector() {}

  //! get the two best electron based on highest pT
  std::pair<int,int> bestByPt();
  //! get the two best electron based on highest SC energy
  std::pair<int,int> bestBySCenergy();
  //! get the two best electron based on best tracker isolation
  std::pair<int,int> bestByTrackerIsolation();
  //! get the two best electron based on best ecal isolation
  std::pair<int,int> bestByEcalIsolation();
  //! get the two best electron based on best hcal isolation
  std::pair<int,int> bestByHcalIsolation();
  //! get the two best electron based on best electron ID LH
  std::pair<int,int> bestByElectronIdLH();

protected:

  std::vector<ElectronQualityData> m_electrons;

  //! compare the pT's
  struct PtSorting : public std::binary_function<ElectronQualityData, ElectronQualityData, bool> {
    bool operator() (ElectronQualityData ele1, ElectronQualityData ele2) 
    { return (ele1.ptAtVtx > ele2.ptAtVtx); }
  };

  //! compare the SC energy
  struct SCenergySorting : public std::binary_function<ElectronQualityData, ElectronQualityData, bool> {
    bool operator() (ElectronQualityData ele1, ElectronQualityData ele2) 
    { return (ele1.SCenergy > ele2.SCenergy); }
  };

  //! compare the tracker isolations
  struct TkIsolationSorting : public std::binary_function<ElectronQualityData, ElectronQualityData, bool> {
    bool operator() (ElectronQualityData ele1, ElectronQualityData ele2) 
    { return (ele1.trackerSumPt < ele2.trackerSumPt); }
  };

  //! compare the ecal isolations
  struct EcalIsolationSorting : public std::binary_function<ElectronQualityData, ElectronQualityData, bool> {
    bool operator() (ElectronQualityData ele1, ElectronQualityData ele2) 
    { return (ele1.ecalSumEt < ele2.ecalSumEt); }
  };

  //! compare the hcal isolations
  struct HcalIsolationSorting : public std::binary_function<ElectronQualityData, ElectronQualityData, bool> {
    bool operator() (ElectronQualityData ele1, ElectronQualityData ele2) 
    { return (ele1.hcalSumEt < ele2.hcalSumEt); }
  };

  //! compare the electron LH
  struct LHSorting : public std::binary_function<ElectronQualityData, ElectronQualityData, bool> {
    bool operator() (ElectronQualityData ele1, ElectronQualityData ele2) 
    { return (ele1.electronIdLH > ele2.electronIdLH); }
  };


};

#endif
