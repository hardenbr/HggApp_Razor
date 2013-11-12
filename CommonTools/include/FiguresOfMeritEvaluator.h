#ifndef FiguresOfMeritEvaluator_h
#define FiguresOfMeritEvaluator_h

#include <vector>
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"

class FiguresOfMeritEvaluator {

 public:
  //! ctor and dtor
  FiguresOfMeritEvaluator();
  ~FiguresOfMeritEvaluator() {}
  
  //! add signal
  void addSignal(const char *nameVar, TH1F* sig);
  void addSignal(const char *nameVar, TH2F* sig);

  //! add backgrounds, at least 1
  void addBackgrounds(TH1F* bkg0, TH1F* bkg1=0, TH1F* bkg2=0, TH1F* bkg3=0);
  void addBackgrounds(TH2F* bkg0, TH2F* bkg1=0, TH2F* bkg2=0, TH2F* bkg3=0);

  //! set the direction of the cut (var>x0 or var<x0)
  void setCutDirection(const char *dir="") { m_direction.push_back(dir); }
  void setCutDirectionXY(const char *dir="") { m_directionXY.push_back(dir); }

  //! draw the comparison for all the variables
  //! option = 0 => sig eff vs. bkg rej
  //! option = 1 => sig eff vs. bkg eff
  void drawResults(const char *fileName, int option=0);

  //! set the title of the signal (default is "signal")
  void setSignalTitle(const char *title="signal") { m_signalTitle = title; }

  //! set the title of the background (default is "background")
  void setBackgroundTitle(const char *title="background") { m_backgroundTitle = title; }

  //! zoom the X-Y axis
  void setRange(double xmin, double xmax, double ymin, double ymax);

 private:

  //! compute FOM vor a certain variable
  TGraphErrors* getFOM1D(const char *nameVar, int option=0);
  TGraphErrors* getFOM2D(const char *nameVar, int option=0);

  //! put a nice axis font
  void AxisFonts(TAxis*  axis, TString coordinate, TString title);

  //! list of variables
  std::vector<TH1F*> m_signalHisto;
  std::vector<TH1F*> m_bkgHisto;
  std::vector<TString> m_names;

  //! list of 2D variables
  std::vector<TH2F*> m_signalHisto2D;
  std::vector<TH2F*> m_bkgHisto2D;
  std::vector<TString> m_names2D;  

  //! direction of the cut (1D)
  std::vector<const char*> m_direction;

  //! direction of the cut (2D)
  std::vector<const char*> m_directionXY;

  //! titles for legend
  const char *m_signalTitle;
  const char *m_backgroundTitle;

  //! ranges of the final plot
  double m_xmin, m_xmax, m_ymin, m_ymax; 

  //! cosmetics
  Int_t   _labelFont      ;
  Float_t _labelOffset    ;
  Float_t _axisLabelSize  ;
  Float_t _titleOffset    ;

};

#endif
