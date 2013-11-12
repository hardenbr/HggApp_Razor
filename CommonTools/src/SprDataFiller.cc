#include <fstream>

#include "CommonTools/include/SprDataFiller.hh"

SprDataFiller::SprDataFiller() {
  _name = "sprtraindata";
  _initialized = false;
  _nentries = -1;
}



void SprDataFiller::setName(const char *name) {
  _name = std::string(name);
}



void SprDataFiller::add(std::string varName, float *var) {

  _data.push_back( std::make_pair( varName, var ) );

}



void SprDataFiller::fill(int ievent, int signal) {

  if( !_initialized ) {
    initialize();
    fillEvent(ievent, signal);
  }
  else {
    fillEvent(ievent, signal);
  }

}



void SprDataFiller::initialize() {

  std::string trainName = _name+"-train.dat";
  std::string validName = _name+"-valid.dat";


  _trainFile.open( trainName.c_str() );
  _validFile.open( validName.c_str() );
  
  int nVariables = _data.size();

  _trainFile << "# SPR training file made by HiggsApp" << std::endl;
  _trainFile << "# number of entries total = " << _nentries << std::endl;
  _trainFile << nVariables << std::endl;

  _validFile << "# SPR validation file made by HiggsApp" << std::endl;
  _validFile << "# number of entries total = " << _nentries << std::endl;
  _validFile << nVariables << std::endl;
  
  std::vector< std::pair< std::string, float* > >::const_iterator dataItr;

  for ( dataItr=_data.begin(); dataItr!=_data.end(); ++dataItr ) {
    _trainFile << dataItr->first << "\t";
    _validFile << dataItr->first << "\t";
  }

  _trainFile << std::endl;
  _validFile << std::endl;
  
  _initialized = true;

}



void SprDataFiller::fillEvent(int ievent, int signal) {
  
  
  if( ((float)ievent)/((float)_nentries) < 0.5 ) {
    
    //    _trainFile << "#\t" << ievent << std::endl;

    std::vector< std::pair < std::string, float* > >::const_iterator dataItr;
    for ( dataItr=_data.begin(); dataItr!=_data.end(); ++dataItr ) {

      _trainFile << *(dataItr->second) << "\t";

    }
    _trainFile << signal << std::endl;
  }
  else {
    
    //    _validFile << "#\t" << ievent << std::endl;

    std::vector< std::pair < std::string, float* > >::const_iterator dataItr;
    for ( dataItr=_data.begin(); dataItr!=_data.end(); ++dataItr ) {

      _validFile << *(dataItr->second) << "\t";

    }
    _validFile << signal << std::endl;

  }

}
