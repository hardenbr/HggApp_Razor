//-------------------------------------------------------
// Description:
//    Class to fill dataset for StatPatternRecognition 
//    classifier training
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef SPRDATAFILLER_H
#define SPRDATAFILLER_H

#include <string>
#include <vector>
#include <fstream>

class SprDataFiller {

private:

  int _nentries;
  std::vector< std::pair < std::string, float* > > _data;
  std::string _name; 
  bool _initialized;
  std::ofstream _trainFile, _validFile;

  void initialize();
  void fillEvent(int ievent, int signal);

public:

  SprDataFiller();
  virtual ~SprDataFiller() {};

  //! set the number of entries of the input tree
  void setEntries(int nentries) { _nentries = nentries; }
  //! set name of the output
  void setName(const char *name);
  //! add a variable
  void add(std::string varName, float *var);
  //! signal = 1; background = 0
  void fill(int ievent, int signal);
  
};

#endif
