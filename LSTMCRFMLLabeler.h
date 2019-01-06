/*
 * Labeler.h
 *
 *  Created on: Mar 16, 2015
 *      Author: mszhang
 */

#ifndef SRC_NNCRF_H_
#define SRC_NNCRF_H_


#include "LibN3L-master/N3L.h"
#include "basic/LSTMCRFMLClassifier.h"
#include "Options.h"
#include "Instance.h"
#include "Example.h"
#include <time.h>


#include "Pipe.h"
#include "Utf.h"

using namespace nr;
using namespace std;

class Labeler {

public:
  std::string nullkey;
  std::string unknownkey;
  std::string seperateKey;

public:
  Alphabet m_featAlphabet;
  Alphabet m_labelAlphabet;

  Alphabet m_wordAlphabet;

  NRVec<Alphabet> m_tagAlphabets;

public:
  Options m_options;

  Pipe m_pipe;

#if USE_CUDA==1
  LSTMCRFMLClassifier<gpu> m_classifier;
#else
  LSTMCRFMLClassifier<cpu> m_classifier;
#endif

public:
  Labeler();
  virtual ~Labeler();

public:

  int createAlphabet(const vector<string>& vecTrainFiles, const vector<vector<size_t> >& vecTrainInstsPositions, const int& tagNum);

  void extractLinearFeatures(vector<string>& features, const Instance* pInstance, int idx);
  void extractFeature(Feature& feat, const Instance* pInstance, int idx);

  void convert2Example(const Instance* pInstance, Example& exam);
  void initialExamples(const vector<Instance>& vecInsts, vector<Example>& vecExams);

public:
  void train(const string& trainFiles, const string& trainCorpusSize, const string& devFiles, const string& testFiles, const string& modelFile, const string& optionFile);
  int predict(const vector<Feature>& features, vector<string>& outputs, const vector<string>& words);
  void test(const string& testFile, const string& outputFile, const string& modelFile);

  void writeModelFile(const string& outputModelFile);
  void loadModelFile(const string& inputModelFile);

};

#endif /* SRC_NNCRF_H_ */
