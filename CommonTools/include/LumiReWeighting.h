#ifndef LumiReWeighting_h
#define LumiReWeighting_h

#include "TH1.h"
#include "TFile.h"
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

class LumiReWeighting {
 public:

  LumiReWeighting( std::string generatedFile, std::string dataFile, std::string GenHistName, std::string DataHistName);
  
  LumiReWeighting( std::vector< float > MC_distr, std::vector< float > Lumi_distr);
  
  void weightOOT_init();
  
  double ITweight( int npv );

  double ITweight3BX( float ave_int );

  double weightOOT( int npv_in_time, int npv_m50nsBX );

  protected:

    std::string generatedFileName_;
    std::string dataFileName_;
    std::string GenHistName_;
    std::string DataHistName_;
    TFile *generatedFile_;
    TFile *dataFile_;
    TH1  *weights_;

    double WeightOOTPU_[25][25];

    bool FirstWarning_;


};


#endif
