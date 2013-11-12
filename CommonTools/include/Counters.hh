//-------------------------------------------------------
// Description:
//    Class for evaluation and display of efficiencies
// Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
//-------------------------------------------------------

#ifndef COUNTERS_H
#define COUNTERS_H

#include <vector>
#include <string>
using namespace std;
class Counters {
 private:
  vector<string> _names;
  vector<float> _counts;
  int _nevents;
  string _title;
 public:
  Counters() {};
  Counters(string title) {_title=title;};
  void SetTitle(const char* title);
  void AddVar(string name);
  void IncrVar(string name);
  void IncrVar(string name, float weight);
  void IncrNevents();
  void Draw();
  void Draw(string name1,string name2);
  float GetVar(string name);
  void Save(const char* filename, const char* option);
  
};




#endif
