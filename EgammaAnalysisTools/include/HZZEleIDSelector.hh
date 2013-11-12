#ifndef HZZEleIDSelector_H
#define HZZEleIDSelector_H

class HZZEleIDSelector {
public:

  HZZEleIDSelector() {}
  ~HZZEleIDSelector() {}
  enum wpfulliso {
    kWP95 = 0,
    kWP90, kWP85, kWP80, kWP70, kWPHWW, kWPHZZ,
    kMVALoose, kMVATight
  };

  enum wpchiso {
    kWP95ChIso = 0,
    kWP90ChIso, kWP85ChIso, kWP80ChIso, kWP70ChIso
  };

  enum mvatype {
    kMVABiased = 0,
    kMVAUnbiased
  };

  enum cutblock {
    all = 0,
    idonly,
    isoonly
  };

  bool output(float pt, float eta, float bdt, float iso, 
	      HZZEleIDSelector::wpfulliso WP, HZZEleIDSelector::mvatype type, HZZEleIDSelector::cutblock cuts = all);
  bool output(float pt, float eta, float bdt, float chiso, 
	      HZZEleIDSelector::wpchiso WP, HZZEleIDSelector::mvatype type, HZZEleIDSelector::cutblock cuts = all);

private:
  int etabin(float eta);
  int ptbinTrg(float pt);
  int ptbinNoTrg(float pt);

};

#endif
