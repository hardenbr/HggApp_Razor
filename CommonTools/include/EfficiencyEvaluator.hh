#ifndef EFFICIENCYEVALUATOR_H
#define EFFICIENCYEVALUATOR_H

#include <vector>
#include <TH1F.h>
#include <TFile.h>

class EfficiencyEvaluator {
private:
  std::vector<TH1F*> _numerators;
  TH1F *_denominator;
  std::vector<TH1F*> _efficiencies;
  TFile *_file;
  const char *m_namefile;
  const char *m_xtitle;
  const char *m_ytitle;
  const char *m_title;
  float m_ymin;
  float m_ymax;

public:
  EfficiencyEvaluator(const char* namefile, const char* opt="RECREATE");
  ~EfficiencyEvaluator();
  //! add a histogram as efficiency numerator
  void AddNumerator(TH1F *numerator) {_numerators.push_back(numerator);};
  //! set the denominator for the finalk efficiency (typically generator)
  void SetDenominator(TH1F *denominator) {_denominator=denominator;};
  //! compute eff wrt previous cut and absolute efficiency
  //! if doAlsoPartialEfficiency=true then compute also efficiency step by step
  void ComputeEfficiencies(bool doAlsoPartialEfficiency=true);
  //! write the efficiency histograms in a ROOT file
  void Write();
  //! draw the cut / total efficiency histograms in a png file
  void DrawSequential();
  //! draw the single cut efficiency (wrt previous cut)
  void DrawSingle();
  //! draw both single and sequential
  void DrawAll();
  //! retrieve the vector of cumulative efficiencies
  std::vector<TH1F*> GetCumulativeEfficiencies();
  //! retrieve the vector of weighted average of cumulative efficiencies
  //! if called with default arguments it integrates over the whole range
  std::vector<float> GetCumulativeEfficiencyAverages(int firstBin=-1, int lastBin=-1);
  //! retrieve the vector of single cut efficiencies
  std::vector<TH1F*> GetSingleCutEfficiencies();
  //! set the X axis title
  void SetXaxisTitle(const char *title) { m_xtitle = title; }
  //! set the Y axis title
  void SetYaxisTitle(const char *title) { m_ytitle = title; }
  //! set the histogram title
  void SetTitle(const char *title) { m_title = title; }
  //! set the histogram minimum
  void SetYaxisMin(float min) { m_ymin = min; }
  //! set the histogram maximum
  void SetYaxisMax(float max) { m_ymax = max; }

};
#endif
