// note: this is the 2012 version from:
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification
// only the ID part is implemented (not the iso and conv rejection)
#ifndef eIDSimpleCutsSelector_H
#define eIDSimpleCutsSelector_H

class eIDSimpleCutsSelector {
public:
  
  eIDSimpleCutsSelector(float deta, float dphi, float see, float hoe, 
			float d0, float dz, float IoEmIoP);
  ~eIDSimpleCutsSelector() {}

  enum wp {
    kVeto = 0,
    kLoose,
    kMedium,
    kTight
  };

  bool output(float eta, int WP);

private:
  float deta_, dphi_, see_, hoe_, d0_, dz_, IoEmIoP_;
};

#endif

