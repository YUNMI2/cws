#ifndef _PARSER_OPTIONS_
#define _PARSER_OPTIONS_

#pragma once

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include "LibN3L-master/N3L.h"

using namespace std;

class Options {
public:
  
  int trainSampleSize;
  int othertrainSampleSize;

  int wordCutOff;
  int featCutOff;
  int spellCutOff;
  int charCutOff;
  int tagCutOff;
  dtype initRange;
  int maxIter;
  int stopIter;
  int batchSize;
  dtype adaEps;
  dtype adaAlpha;
  dtype regParameter;
  dtype dropProb;
  dtype belta1;
  dtype belta2;
  

  int word_M_dim;
  int word_N_dim;
  string kernel_size;
  string stride;
  int kernel_size_n;
  int kernel_size_m;
  int stride_n;
  int stride_m;
  int cnnHiddenSize;
  int cnnOutputSize;
  int n_cnn_filters;
  int linearHiddenSize;
  int hiddenSize;
  int rnnHiddenSize;
  int wordEmbSize;
  int wordcontext;
  bool wordEmbFineTune;
  int tagEmbSize;
  bool tagEmbFineTune;
	
  int spellEmbSize;
  bool spellEmbFineTune;

  int charEmbSize;
  int charcontext;
  bool charEmbFineTune;
  int charhiddenSize;

  int verboseIter;
  bool saveIntermediate;
  bool train;
  int maxInstance;
  vector<string> testFiles;
  string outBest;
  bool seg;
  int relu;
  int atomLayers;
  int rnnLayers;

  int spellrnnHiddenSize;
  int spellrnnOutputsize;

  Options() {
    trainSampleSize = 1000;  
	othertrainSampleSize = 1000;

    wordCutOff = 0;
    featCutOff = 0;
    spellCutOff = 0;
    charCutOff = 0;
    tagCutOff = 0;
    initRange = 0.01;
    maxIter = 1000;
    stopIter = 50;
    batchSize = 1;
    adaEps = 1e-6;
    adaAlpha = 0.01;
    regParameter = 1e-8;
    dropProb = 0.0;
    belta1 = 0.9;
    belta2 = 0.999;

  
    word_M_dim = 0;
    word_N_dim = 0;

    cnnHiddenSize = 0;
    cnnOutputSize = 0;
    n_cnn_filters = 0;
    kernel_size_n = 0;
    kernel_size_m = 0;
    stride_n = 0;
    stride_m = 0;




    linearHiddenSize = 30;
    hiddenSize = 200;
    rnnHiddenSize = 300;
    wordEmbSize = 50;
    wordcontext = 2;
    wordEmbFineTune = true;
    tagEmbSize = 50;
    tagEmbFineTune = true;
    
    spellEmbSize = 50;
    spellEmbFineTune = true;
    
    charEmbSize = 50;
    charcontext = 2;
    charEmbFineTune = true;
    charhiddenSize = 50;

    verboseIter = 100;
    saveIntermediate = true;
    train = false;
    maxInstance = -1;
    testFiles.clear();
    outBest = "";
    relu = 0;
    seg = false;
    atomLayers = 1;
    rnnLayers = 1;
    spellrnnHiddenSize = 50;
    spellrnnOutputsize = 50;

  }

  virtual ~Options() {

  }

  void setOptions(const vector<string> &vecOption) {
    int i = 0;
    for (; i < vecOption.size(); ++i) {
      pair<string, string> pr;
      string2pair(vecOption[i], pr, '=');
	  if (pr.first == "trainSampleSize")
		trainSampleSize = atoi(pr.second.c_str());
	  if (pr.first == "othertrainSampleSize")
		othertrainSampleSize = atoi(pr.second.c_str());
      
	  if (pr.first == "wordCutOff")
        wordCutOff = atoi(pr.second.c_str());
      if (pr.first == "featCutOff")
        featCutOff = atoi(pr.second.c_str());
      if (pr.first == "spellCutOff")
        spellCutOff = atoi(pr.second.c_str());
      if (pr.first == "charCutOff")
        charCutOff = atoi(pr.second.c_str());
      if (pr.first == "tagCutOff")
        tagCutOff = atoi(pr.second.c_str());        
      if (pr.first == "initRange")
        initRange = atof(pr.second.c_str());
      if (pr.first == "maxIter")
        maxIter = atoi(pr.second.c_str());
      if (pr.first == "stopIter")
        stopIter = atoi(pr.second.c_str());
      if (pr.first == "batchSize")
        batchSize = atoi(pr.second.c_str());
      if (pr.first == "adaEps")
        adaEps = atof(pr.second.c_str());
      if (pr.first == "adaAlpha")
        adaAlpha = atof(pr.second.c_str());
      if (pr.first == "regParameter")
        regParameter = atof(pr.second.c_str());
      if (pr.first == "dropProb")
        dropProb = atof(pr.second.c_str());
      if (pr.first == "belta1")
        belta1 = atof(pr.second.c_str());
      if (pr.first == "belta2")
        belta2 = atof(pr.second.c_str());

      if (pr.first == "linearHiddenSize")
        linearHiddenSize = atoi(pr.second.c_str());
      if (pr.first == "hiddenSize")
        hiddenSize = atoi(pr.second.c_str());
      if (pr.first == "rnnHiddenSize")
        rnnHiddenSize = atoi(pr.second.c_str());
      if (pr.first == "wordcontext")
        wordcontext = atoi(pr.second.c_str());
      if (pr.first == "wordEmbSize")
        wordEmbSize = atoi(pr.second.c_str());
      if (pr.first == "wordEmbFineTune")
        wordEmbFineTune = (pr.second == "true") ? true : false;
      if (pr.first == "tagEmbSize")
        tagEmbSize = atoi(pr.second.c_str());
      if (pr.first == "tagEmbFineTune")
        tagEmbFineTune = (pr.second == "true") ? true : false;        	
      if (pr.first == "spellEmbSize")
        spellEmbSize = atoi(pr.second.c_str());
      if (pr.first == "spellEmbFineTune")
        spellEmbFineTune = (pr.second == "true") ? true : false;        	
      if (pr.first == "charcontext")
        charcontext = atoi(pr.second.c_str());
      if (pr.first == "charEmbSize")
        charEmbSize = atoi(pr.second.c_str());
      if (pr.first == "charEmbFineTune")
        charEmbFineTune = (pr.second == "true") ? true : false;
      if (pr.first == "charhiddenSize")
        charhiddenSize = atoi(pr.second.c_str());
        
      if (pr.first == "verboseIter")
        verboseIter = atoi(pr.second.c_str());
      if (pr.first == "train")
        train = (pr.second == "true") ? true : false;
      if (pr.first == "saveIntermediate")
        saveIntermediate = (pr.second == "true") ? true : false;
      if (pr.first == "maxInstance")
        maxInstance = atoi(pr.second.c_str());
      if (pr.first == "testFile")
        testFiles.push_back(pr.second);
      if (pr.first == "outBest")
        outBest = pr.second;
      if (pr.first == "relu")
        relu = atoi(pr.second.c_str());
      if (pr.first == "seg")
        seg = (pr.second == "true") ? true : false;
      if (pr.first == "atomLayers")
        atomLayers = atoi(pr.second.c_str());
      if (pr.first == "rnnLayers")
        rnnLayers = atoi(pr.second.c_str());
      if (pr.first == "kernel_size_n")
        kernel_size_n = atoi(pr.second.c_str());
      if (pr.first == "kernel_size_m")
        kernel_size_m = atoi(pr.second.c_str());
      if (pr.first == "stride_n")
        stride_n = atoi(pr.second.c_str());
      if (pr.first == "stride_m")
        stride_m = atoi(pr.second.c_str());
      if (pr.first == "word_M_dim")
        word_M_dim = atoi(pr.second.c_str());
      if (pr.first == "word_N_dim")
        word_N_dim = atoi(pr.second.c_str());
      if (pr.first == "cnnHiddenSize")
        cnnHiddenSize = atoi(pr.second.c_str());
      if (pr.first == "cnnOutputSize")
        cnnOutputSize = atoi(pr.second.c_str());
      if (pr.first == "n_cnn_filters")
        n_cnn_filters = atoi(pr.second.c_str());
      if (pr.first == "spellrnnHiddenSize")
        spellrnnHiddenSize = atoi(pr.second.c_str());
      if (pr.first == "spellrnnOutputsize")
        spellrnnOutputsize = atoi(pr.second.c_str());
     //if (pr.first == "stride_m")
     //tride_m = atoi(pr.second.c_str());

    }
  }

  void showOptions() {
    std::cout << "trainSampleSize = " << trainSampleSize << std::endl;
    std::cout << "othertrainSampleSize = " << othertrainSampleSize << std::endl;
	
    std::cout << "wordCutOff = " << wordCutOff << std::endl;
    std::cout << "featCutOff = " << featCutOff << std::endl;
    std::cout << "spellCutOff = " << spellCutOff << std::endl;
    std::cout << "charCutOff = " << charCutOff << std::endl;
    std::cout << "tagCutOff = " << tagCutOff << std::endl;
    std::cout << "initRange = " << initRange << std::endl;
    std::cout << "maxIter = " << maxIter << std::endl;
    std::cout << "stopIter = " << stopIter << std::endl;
    std::cout << "batchSize = " << batchSize << std::endl;
    std::cout << "adaEps = " << adaEps << std::endl;
    std::cout << "adaAlpha = " << adaAlpha << std::endl;
    std::cout << "regParameter = " << regParameter << std::endl;
    std::cout << "dropProb = " << dropProb << std::endl;
    std::cout << "belta1 = " << belta1 << std::endl;
    std::cout << "belta2 = " << belta2 << std::endl;


    std::cout << "word_M_dim = " << word_M_dim << std::endl;
    std::cout << "word_N_dim = " << word_N_dim << std::endl;
    std::cout << "kernel_size_n = " << kernel_size_n << std::endl;
    std::cout << "kernel_size_m = " << kernel_size_m << std::endl;
    std::cout << "stride_n = " << stride_n << std::endl;
    std::cout << "stride_m = " << stride_m << std::endl;
    std::cout << "cnnHiddenSize = " <<  cnnHiddenSize << std::endl;
    std::cout << "cnnOutputSize = " << cnnOutputSize << std::endl;
    std::cout << "n_cnn_filters = " << n_cnn_filters << std::endl;


    std::cout << "linearHiddenSize = " << linearHiddenSize << std::endl;
    std::cout << "hiddenSize = " << hiddenSize << std::endl;
    std::cout << "rnnHiddenSize = " << rnnHiddenSize << std::endl;
    std::cout << "wordEmbSize = " << wordEmbSize << std::endl;
    std::cout << "wordcontext = " << wordcontext << std::endl;
    std::cout << "wordEmbFineTune = " << wordEmbFineTune << std::endl;
    std::cout << "tagEmbSize = " << tagEmbSize << std::endl;
    std::cout << "tagEmbFineTune = " << tagEmbFineTune << std::endl;
    std::cout << "spellEmbSize = " << spellEmbSize << std::endl;
    std::cout << "spellEmbFineTune = " << spellEmbFineTune << std::endl;
    std::cout << "charEmbSize = " << charEmbSize << std::endl;
    std::cout << "charcontext = " << charcontext << std::endl;
    std::cout << "charEmbFineTune = " << charEmbFineTune << std::endl;
    std::cout << "charhiddenSize = " << charhiddenSize << std::endl;

    std::cout << "verboseIter = " << verboseIter << std::endl;
    std::cout << "saveItermediate = " << saveIntermediate << std::endl;
    std::cout << "train = " << train << std::endl;
    std::cout << "maxInstance = " << maxInstance << std::endl;
    for (int idx = 0; idx < testFiles.size(); idx++) {
      std::cout << "testFile = " << testFiles[idx] << std::endl;
    }
    std::cout << "outBest = " << outBest << std::endl;
    std::cout << "relu = " << relu << std::endl;
    std::cout << "seg = " << seg << std::endl;
    std::cout << "atomLayers = " << atomLayers << std::endl;
    std::cout << "rnnLayers = " << rnnLayers << std::endl;
    std::cout << "spellrnnHiddenSize = " << spellrnnHiddenSize << std::endl;
    std::cout << "spellrnnOutputsize = " << spellrnnOutputsize << std::endl;
  }

  void load(const std::string& infile) {
    ifstream inf;
    inf.open(infile.c_str());
    vector<string> vecLine;
    while (1) {
      string strLine;
      if (!my_getline(inf, strLine)) {
        break;
      }
      if (strLine.empty())
        continue;
      vecLine.push_back(strLine);
    }
    inf.close();
    setOptions(vecLine);
  }

  void writeModel(LStream &outf) {


    WriteVector(outf, testFiles);
    WriteString(outf, outBest);
	
    WriteBinary(outf, trainSampleSize);
    WriteBinary(outf, othertrainSampleSize);

    WriteBinary(outf, wordCutOff);
    WriteBinary(outf, featCutOff);
    WriteBinary(outf, spellCutOff);
    WriteBinary(outf, charCutOff);
    WriteBinary(outf, tagCutOff);
    WriteBinary(outf, initRange);
    WriteBinary(outf, maxIter);
    WriteBinary(outf, stopIter);
    WriteBinary(outf, batchSize);
    WriteBinary(outf, adaEps);
    WriteBinary(outf, adaAlpha);
    WriteBinary(outf, regParameter);
    WriteBinary(outf, dropProb);
    WriteBinary(outf, belta1);
    WriteBinary(outf, belta2);
    WriteBinary(outf, linearHiddenSize);
    WriteBinary(outf, hiddenSize);
    WriteBinary(outf, rnnHiddenSize);
    WriteBinary(outf, wordEmbSize);
    WriteBinary(outf, wordcontext);
    WriteBinary(outf, wordEmbFineTune);
    WriteBinary(outf, tagEmbSize);
    WriteBinary(outf, tagEmbFineTune);
    WriteBinary(outf, spellEmbSize);
    WriteBinary(outf, spellEmbFineTune);
    WriteBinary(outf, charEmbSize);
    WriteBinary(outf, charcontext);
    WriteBinary(outf, charEmbFineTune);
    WriteBinary(outf, charhiddenSize);
    WriteBinary(outf, verboseIter);
    WriteBinary(outf, saveIntermediate);
    WriteBinary(outf, train);
    WriteBinary(outf, maxInstance);
    WriteBinary(outf, seg);
    WriteBinary(outf, relu);
    WriteBinary(outf, atomLayers);
    WriteBinary(outf, rnnLayers);
    WriteBinary(outf, cnnHiddenSize);
    WriteBinary(outf, cnnOutputSize);
    WriteBinary(outf, n_cnn_filters);
    WriteBinary(outf, kernel_size_n);
    WriteBinary(outf, kernel_size_m);
    WriteBinary(outf, stride_n);
    WriteBinary(outf, stride_m);
    WriteBinary(outf, word_M_dim);
    WriteBinary(outf, word_N_dim);
    WriteBinary(outf, spellrnnHiddenSize);
    WriteBinary(outf, spellrnnOutputsize);
 //   WriteBinary(outf, rnnLayers);

  }

  void loadModel(LStream &inf) {
    ReadVector(inf, testFiles);
    ReadString(inf, outBest);



   /* ReadBinary(inf,cnnHiddenSize);
    ReadBinary(inf,cnnOutputSize);
    ReadBinary(inf,n_cnn_filters);
    ReadBinary(inf,kernel_size_n);
    ReadBinary(inf,kernel_size_m);
    ReadBinary(inf,stride_n);
    ReadBinary(inf,stride_m);*/
   // ReadBinary(inf,cnnHiddenSize);
    ReadBinary(inf, trainSampleSize);
    ReadBinary(inf, othertrainSampleSize);

    ReadBinary(inf, wordCutOff);
    ReadBinary(inf, featCutOff);
    ReadBinary(inf, spellCutOff);
    ReadBinary(inf, charCutOff);
    ReadBinary(inf, tagCutOff);
    ReadBinary(inf, initRange);
    ReadBinary(inf, maxIter);
    ReadBinary(inf, stopIter);
    ReadBinary(inf, batchSize);
    ReadBinary(inf, adaEps);
    ReadBinary(inf, adaAlpha);
    ReadBinary(inf, regParameter);
    ReadBinary(inf, dropProb);
    ReadBinary(inf, belta1);
    ReadBinary(inf, belta2);

    ReadBinary(inf, linearHiddenSize);
    ReadBinary(inf, hiddenSize);
    ReadBinary(inf, rnnHiddenSize);
    ReadBinary(inf, wordEmbSize);
    ReadBinary(inf, wordcontext);
    ReadBinary(inf, wordEmbFineTune);
    ReadBinary(inf, tagEmbSize);
    ReadBinary(inf, tagEmbFineTune);
    ReadBinary(inf, spellEmbSize);
    ReadBinary(inf, spellEmbFineTune);
    ReadBinary(inf, charEmbSize);
    ReadBinary(inf, charcontext);
    ReadBinary(inf, charEmbFineTune);
    ReadBinary(inf, charhiddenSize);
    ReadBinary(inf, verboseIter);
    ReadBinary(inf, saveIntermediate);
    ReadBinary(inf, train);
    ReadBinary(inf, maxInstance);
    ReadBinary(inf, seg);
    ReadBinary(inf, relu);
    ReadBinary(inf, atomLayers);
    ReadBinary(inf, rnnLayers);
    ReadBinary(inf, cnnHiddenSize);
    ReadBinary(inf, cnnOutputSize);
    ReadBinary(inf, n_cnn_filters);
    ReadBinary(inf, kernel_size_n);
    ReadBinary(inf, kernel_size_m);
    ReadBinary(inf, stride_n);
    ReadBinary(inf, stride_m);
    ReadBinary(inf, word_M_dim);
    ReadBinary(inf, word_N_dim);
    ReadBinary(inf, spellrnnHiddenSize);
    ReadBinary(inf, spellrnnOutputsize);
  }

  
};

#endif

