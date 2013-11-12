#ifndef SimpleCutsIDSelector_H
#define SimpleCutsIDSelector_H

class SimpleCutsIDSelector {
public:

  SimpleCutsIDSelector() {}
  ~SimpleCutsIDSelector() {}
  enum wp {
    kWP95 = 0,
    kWP90, kWP85, kWP80, kWP70
  };

  bool outputid(float eta, float see, float dphi, float deta, float hoe, 
		SimpleCutsIDSelector::wp WP);
  bool outputiso(float eta, float tk, float ecal, float hcal,
		 SimpleCutsIDSelector::wp WP);
  bool outputconv(float eta, int misshits, float dist, float dcot, 
		  SimpleCutsIDSelector::wp WP);
  
private:
  int etabin(float eta);
};

#endif
