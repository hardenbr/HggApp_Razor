#ifndef SKIMMER_H
#define SKIMMER_H

#include <vector>

class Skimmer {

public:
  Skimmer(const char *skimTxtFile);
  virtual ~Skimmer();

  void readFile();
  bool output(int event);


private:
  const char *m_file;
  std::vector<int> bits;

};

#endif
