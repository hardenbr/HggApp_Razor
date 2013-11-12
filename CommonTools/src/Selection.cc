#include "CommonTools/include/Selection.hh"
#include <iostream>
#include <fstream>

Selection::Selection(std::string fileCuts, std::string fileSwitches) {
  _fileCut=fileCuts;
  _fileSwitch=fileSwitches;
}
Selection::~Selection() {}

std::pair<float,float> Selection::readIntervalFromFile(std::string name) {

  std::pair<float,float> defRange(-999999.,999999.);
  float min, max;
  std::ifstream setfile(_fileCut.c_str());
  if(!setfile.good()) {
    std::cout << "Selection::Error!   Unable to open the file " << _fileCut << std::endl;
    return defRange;
  }
  else {
    std::string var; 
    bool found=false;
    while (1) {
      setfile >> var >> min >> max;
      if(!setfile.good()) break;
      if(var==name) { 
	found=true;
	break;
      }
    }
    if(!found) {
      std::cout<< "Selection::Warning! cannot read cut settings for variable " 
	       << name << " in file " << _fileCut << std::endl;
      return defRange;
    }
  }
  setfile.close();
  std::pair<float,float> range(min,max);
  return range;
}

int Selection::readSwitchFromFile(std::string name) {

  int Switch;
  std::ifstream setfile(_fileSwitch.c_str());
  if(!setfile.good()) {
    std::cout << "Selection::Error!   Unable to open the file " << _fileSwitch << std::endl;
    return 0;
  }
  else {
    std::string var; 
    bool found=false;
    while (1) {
      setfile >> var >> Switch;
      if(!setfile.good()) break;
      if(var==name) { 
	found=true;
	break;
      }
    }
    if(!found) {
      std::cout<< "Selection::Warning! cannot read settings for switch " 
	       << name << " in file " << _fileSwitch
	       << "switching it to false" << std::endl;
      return 0;
    }
  }
  setfile.close();
  return Switch;
}

std::string Selection::readStringFromFile(std::string name) {

  std::string StringPar;
  std::ifstream setfile(_fileSwitch.c_str());
  if(!setfile.good()) {
    std::cout << "Selection::Error!   Unable to open the file " << _fileSwitch << std::endl;
    return 0;
  }
  else {
    std::string var; 
    bool found=false;
    while (1) {
      setfile >> var >> StringPar;
      if(!setfile.good()) break;
      if(var==name) { 
	found=true;
	break;
      }
    }
    if(!found) {
      std::cout<< "Selection::Warning! cannot read settings for string parameter " 
	       << name << " in file " << _fileSwitch
	       << "switching it to empty string" << std::endl;
      return std::string("");
    }
  }
  setfile.close();
  return StringPar;
}


void Selection::addCut(std::string name) {
  std::pair<float,float> range = readIntervalFromFile(name);
  _cut.insert( make_pair( name, range) );

  // add the correspondent switches
  int switchStatus = (range.second < 999998 ) ? readSwitchFromFile(name) : 0;
  _switch.insert( make_pair( name, switchStatus) );
}

void Selection::addSwitch(std::string name) {
  int Switch = readSwitchFromFile(name);
  _switch.insert( make_pair( name, Switch) );
}

void Selection::addStringParameter(std::string name) {
  std::string StringPar = readStringFromFile(name);
  _stringPar.insert( make_pair( name, StringPar) );
}

std::pair<float,float> Selection::getCut(std::string name) {
  return _cut[name];
} 

int Selection::getSwitch(std::string name) {
  return _switch[name];
}

std::string Selection::getStringParameter(std::string name) {
  return _stringPar[name];
}

bool Selection::passCut(std::string name, float var) {
  std::pair<float,float> range = _cut[name];
  return (var>=range.first && var<=range.second);
}

bool Selection::passCut(std::string name, std::vector<float> vars) {
  std::pair<float,float> range = _cut[name];
  std::vector<float>::iterator iter;
  for(iter=vars.begin();iter!=vars.end();iter++){
    if(*iter<=range.first || *iter>=range.second) return false;
  }
  return true;
}

float Selection::getUpperCut(std::string name) {
  std::pair<float,float> range = _cut[name];
  return range.second;
}

float Selection::getLowerCut(std::string name) {
  std::pair<float,float> range = _cut[name];
  return range.first;
}

void Selection::summary() {
  std::map<std::string,std::pair<float,float> >::iterator iter;
  for( iter = _cut.begin(); iter!=_cut.end(); iter++) {
    std::pair<float,float> range = iter->second;
    std::string cutStatus;
    std::map<std::string, int>::iterator switchItr = _switch.find(iter->first);
    // if not found in the switches list, default of the cut is ON
    const char* e0;
    const char* en="\033[0m";
    if(switchItr==_switch.end() || switchItr->second==1) {
      e0="\033[44;37m";
      cutStatus="ON";
    }
    else {
      e0="\033[41;37m";
      cutStatus="OFF";
    }
    std::cout << "*** " << iter->first << " within [" << range.first << " - " << range.second << "]" 
	      << e0 << "\tStatus is: " << cutStatus << en
	      << std::endl; 
  }
  
  std::map< std::string,std::string >::iterator strpar;
  for( strpar = _stringPar.begin(); strpar!=_stringPar.end(); strpar++) {
    const char* e0 = "\033[42;37m";
    const char* en="\033[0m";
    std::cout << "*** " << e0 << "Parameter " << strpar->first
              << " has value " << strpar->second << en << std::endl;
  }
}

void Selection::append(Selection* selectionToAdd) {

  std::map<std::string,std::pair<float,float> > cutsToAdd = selectionToAdd->getSelection();
  std::map<std::string,int> switchesToAdd = selectionToAdd->getSwitches();

  if(cutsToAdd.size() > 0 ) {
    std::map<std::string,std::pair<float,float> >::iterator iterCut;
    std::map<std::string,int>::iterator iterSwitch=switchesToAdd.begin();
    for( iterCut = cutsToAdd.begin(); iterCut != cutsToAdd.end(); iterCut++ ) {
      _cut.insert( *iterCut );
      _switch.insert( *iterSwitch );
      ++iterSwitch;
    }
  }
  else {
    std::cout << "Selection::Warning!  Appending an empty selection to " 
	      << _fileSwitch << " !" << std::endl;
  }

}
