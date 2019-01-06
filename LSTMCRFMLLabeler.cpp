/*
 * Labeler.cpp
 *
 *  Created on: Mar 16, 2015
 *      Author: mszhang
 */


#include "LSTMCRFMLLabeler.h"

#include "Argument_helper.h"

#include <sstream>

Labeler::Labeler() {
  nullkey = "-null-";
  unknownkey = "-unknown-";
  seperateKey = "#";

}

Labeler::~Labeler() {
  m_classifier.release();
}

int Labeler::createAlphabet(const vector<string>& vecTrainFiles, const vector<vector<size_t> >& vecTrainInstsPositions, const int& tagNum) {
  	cout << "Creating Alphabet..." << endl;

  	int numInstance;
  	hash_map<string, int> feature_stat;
  	hash_map<string, int> word_stat;
  	vector<hash_map<string, int> > tag_stat;
  	m_labelAlphabet.clear();

  	tag_stat.resize(tagNum);
  	m_tagAlphabets.resize(tagNum);

	for(int i_file = 0; i_file < vecTrainFiles.size(); i_file++){
		m_pipe.initInputFile(vecTrainFiles[i_file].c_str());
		size_t posi;
  		for (numInstance = 0; numInstance < vecTrainInstsPositions[i_file].size(); numInstance++) {
			posi = vecTrainInstsPositions[i_file][numInstance];
			Instance *pInstance = m_pipe.readInstanceByPosi(posi);
		
    		const vector<string> &words = pInstance->words;
    		const vector<string> &labels = pInstance->labels;
    		const vector<vector<string> > &sparsefeatures = pInstance->sparsefeatures;

			const vector<vector<string> > &tagfeatures = pInstance->tagfeatures;
    		for (int iter_tag = 0; iter_tag < tagfeatures.size(); iter_tag++) {
      			assert(tagNum == tagfeatures[iter_tag].size());
    		}

    		vector<string> features;
    		int curInstSize = labels.size();
    		int labelId;
    		for (int i = 0; i < curInstSize; ++i) {
			vector<string> vecCurLabels;
			split_bychar(labels[i], vecCurLabels,'/');
			for (int id_label = 0;id_label < vecCurLabels.size();id_label++){
				if (vecCurLabels[id_label].compare("A-SEG") == 0 or vecCurLabels[id_label].compare("a-seg") == 0)//如果label是A-SEG，那么就不考虑
					continue;
      				labelId = m_labelAlphabet.from_string(vecCurLabels[id_label]);

			}
      			string curword = normalize_to_lowerwithdigit(words[i]);
      			word_stat[curword]++;
      			for (int j = 0; j < sparsefeatures[i].size(); j++)
        			feature_stat[sparsefeatures[i][j]]++;
      			for (int j = 0; j < tagfeatures[i].size(); j++)
        			tag_stat[j][tagfeatures[i][j]]++;
    		}

    		if (m_options.maxInstance > 0 && numInstance == m_options.maxInstance)
      			break;
  		}
		m_pipe.uninitInputFile();
	}

  	cout << numInstance << " " << endl;
  	cout << "Label num: " << m_labelAlphabet.size() << endl;
  	cout << "Total word num: " << word_stat.size() << endl;
  	cout << "Total feature num: " << feature_stat.size() << endl;
	// tag print information
  	cout << "tag num = " << tagNum << endl;
  	for (int iter_tag = 0; iter_tag < tagNum; iter_tag++) {
    		cout << "Total tag " << iter_tag << " num: " << tag_stat[iter_tag].size() << endl;
  	}
  	m_featAlphabet.clear();
  	m_wordAlphabet.clear();
 	m_wordAlphabet.from_string(nullkey);
  	m_wordAlphabet.from_string(unknownkey);

  	//tag apheabet init
  	for (int i = 0; i < tagNum; i++) {
    		m_tagAlphabets[i].clear();
    		m_tagAlphabets[i].from_string(nullkey);
    		m_tagAlphabets[i].from_string(unknownkey);
  	}

  	hash_map<string, int>::iterator feat_iter;
  	for (feat_iter = feature_stat.begin(); feat_iter != feature_stat.end(); feat_iter++) {
    		if (feat_iter->second > m_options.featCutOff) {
      			m_featAlphabet.from_string(feat_iter->first);
    		}
  	}

  	for (feat_iter = word_stat.begin(); feat_iter != word_stat.end(); feat_iter++) {
   		if (!m_options.wordEmbFineTune || feat_iter->second > m_options.wordCutOff) {
      			m_wordAlphabet.from_string(feat_iter->first);
    		}
  	}

  	cout << "before tag alphabet line 121" << endl;
	// tag cut off, default tagCutOff is zero
  	for (int i = 0; i < tagNum; i++) {
    		for (feat_iter = tag_stat[i].begin(); feat_iter != tag_stat[i].end(); feat_iter++) {
      			if (!m_options.tagEmbFineTune || feat_iter->second > m_options.tagCutOff) {
        			m_tagAlphabets[i].from_string(feat_iter->first);
      			}
    		}
  	}

  	cout << "Remain feature num: " << m_featAlphabet.size() << endl;
  	cout << "Remain words num: " << m_wordAlphabet.size() << endl;
	// tag Remain num print
  	for (int i = 0; i < tagNum; i++) {
    		cout << "Remain tag " << i << " num: " << m_tagAlphabets[i].size() << endl;
  	}

  	m_labelAlphabet.set_fixed_flag(true);
  	m_featAlphabet.set_fixed_flag(true);
  	m_wordAlphabet.set_fixed_flag(true);

	// tag Alphabet fixed  
  	for (int iter_tag = 0; iter_tag < tagNum; iter_tag++) {
    		m_tagAlphabets[iter_tag].set_fixed_flag(true);
  	}
  
	cout << "Creating alphabet ending!!" << endl;
  	return 0;
}



void Labeler::extractFeature(Feature& feat, const Instance* pInstance, int idx) {
  	feat.clear();

  	const vector<string>& words = pInstance->words;
  	int sentsize = words.size();
  	string curWord = idx >= 0 && idx < sentsize ? normalize_to_lowerwithdigit(words[idx]) : nullkey;

  	// word features
  	int unknownId = m_wordAlphabet.from_string(unknownkey);

  	int curWordId = m_wordAlphabet.from_string(curWord);
  	if (curWordId >= 0)
    		feat.words.push_back(curWordId);
  	else
    		feat.words.push_back(unknownId);

  	// tag features
  	const vector<vector<string> > &tagfeatures = pInstance->tagfeatures;
  	int tagNum = tagfeatures[idx].size();
  	for (int i = 0; i < tagNum; i++) {
    		unknownId = m_tagAlphabets[i].from_string(unknownkey);
    		int curTagId = m_tagAlphabets[i].from_string(tagfeatures[idx][i]);
    		if (curTagId >= 0)
      			feat.tags.push_back(curTagId);
    		else
      			feat.tags.push_back(unknownId);
  	}

  	const vector<string>& linear_features = pInstance->sparsefeatures[idx];
  	for (int i = 0; i < linear_features.size(); i++) {
    		int curFeatId = m_featAlphabet.from_string(linear_features[i]);
    		if (curFeatId >= 0)
      			feat.linear_features.push_back(curFeatId);
 	}

}

void Labeler::convert2Example(const Instance* pInstance, Example& exam) {
  	exam.clear();
  	const vector<string> &labels = pInstance->labels;
  	int curInstSize = labels.size();
  	for (int i = 0; i < curInstSize; ++i) {
		int numLabel1s = m_labelAlphabet.size();
    		string orcale = labels[i];
    		vector<int> curlabels;
		if (orcale.compare("a-seg") == 0 || orcale.compare("A-SEG") == 0){
			for (int j = 0; j < numLabel1s; ++j) {
				curlabels.push_back(1);
			}
			
		} 
		else{
			for (int j = 0; j < numLabel1s; ++j)
				curlabels.push_back(0);
			vector<string> vecCurLabels;
			split_bychar(orcale, vecCurLabels,'/');
			for (int id_label = 0; id_label < vecCurLabels.size(); id_label ++){
    				for (int j = 0; j < numLabel1s; ++j) {
      					string str = m_labelAlphabet.from_id(j);
      					if (str.compare(vecCurLabels[id_label]) == 0)
        					curlabels[j] = 1;
    				}
			}
		}
/*		cout << "labels:" << orcale << endl;
		cout << "curlabels " ;
		for (int tmp_ix = 0; tmp_ix < curlabels.size();tmp_ix ++){
			cout << curlabels[tmp_ix] << " " ;
		}
		cout << endl;*/
    		exam.m_labels.push_back(curlabels);
    		Feature feat;
    		extractFeature(feat, pInstance, i);
    		exam.m_features.push_back(feat);   
  	}
//	exit(0);
}

void Labeler::initialExamples(const vector<Instance>& vecInsts, vector<Example>& vecExams) {
  	int numInstance;
  	for (numInstance = 0; numInstance < vecInsts.size(); numInstance++) {
    		const Instance *pInstance = &vecInsts[numInstance];
    		Example curExam;
    		convert2Example(pInstance, curExam);
    		vecExams.push_back(curExam);

    		if (m_options.maxInstance > 0 && numInstance == m_options.maxInstance)
      			break;
  	}

}


void Labeler::train(const string& trainFiles, const string& trainCorpusSize, const string& devFiles, const string& testFiles, const string& modelFile, const string& optionFile) {
	/*
		xxxFiles:多个xxx语料，如：pku,msr,ctb
		trainCorpusSize:每次迭代中每个训练语料对应的size
		modelFile:模型文件
		optionFile:配置文件
	*/

  	//load option file
	if (optionFile != "")
    	m_options.load(optionFile);
  	m_options.showOptions();

    //get train&&dev&&test filename
	vector<string> vecTrainFiles, vecDevFiles, vecTestFiles, vecTrainCorpusSize;
	split_bychar(trainFiles, vecTrainFiles, ',');
	split_bychar(trainCorpusSize, vecTrainCorpusSize, ',');
	assert(vecTrainFiles.size() == vecTrainCorpusSize.size());

	split_bychar(devFiles, vecDevFiles, ',');
	split_bychar(testFiles, vecTestFiles, ',');
 	assert(vecDevFiles.size() == vecTestFiles.size());
    cout << vecTrainFiles.size() << " " << vecTrainCorpusSize.size() << endl;
    cout << vecDevFiles.size() << " " << vecTestFiles.size() << endl;
  	
    //read file 
    //train: big, read instance position
    //dev&&test: small, read instance
	vector<vector<size_t> > vecTrainInstsPositions(vecTrainFiles.size());
	vector<vector<Instance> > vecDevInsts(vecDevFiles.size()), vecTestInsts(vecTestFiles.size());
  	static vector<Instance> decodeInstResults;
  	static Instance curDecodeInst;
	
	vector<bool> bCurIterBetters;
	vector<dtype> bestDISs;
	for(int i_dev = 0; i_dev < vecDevFiles.size(); i_dev++){
		bCurIterBetters.push_back(false);
		bestDISs.push_back(0);
	}
    cout << "........................Start Reading Instance............................" << endl;
    //load instance
	for(int i = 0; i < vecTrainFiles.size(); i++)
		m_pipe.readInstancesPosi(vecTrainFiles[i], vecTrainInstsPositions[i], m_options.maxInstance);	
	for(int i = 0; i < vecDevFiles.size(); i++)
		m_pipe.readInstances(vecDevFiles[i], vecDevInsts[i], m_options.maxInstance);	
	for(int i = 0; i < vecTestFiles.size(); i++)
		m_pipe.readInstances(vecTestFiles[i], vecTestInsts[i], m_options.maxInstance);	
    
    	
	int tagNum = vecDevInsts[0][0].tagfeatures[0].size();
    cout << "tagNum:" << tagNum << endl;
	createAlphabet(vecTrainFiles, vecTrainInstsPositions, tagNum);

  	NRMat<dtype> wordEmb;
  	wordEmb.resize(m_wordAlphabet.size(), m_options.wordEmbSize);
  	wordEmb.randu(1000);

  	NRVec<NRMat<dtype> > tagEmbs(m_tagAlphabets.size());
  	for (int idx = 0; idx < tagEmbs.size(); idx++){
    	tagEmbs[idx].resize(m_tagAlphabets[idx].size(), m_options.tagEmbSize);
    	tagEmbs[idx].randu(1002 + idx);
  	}
  
	// set parameter
  	m_classifier.setWordEmbFinetune(m_options.wordEmbFineTune);
  	m_classifier.init(wordEmb, m_options.wordcontext, tagEmbs, m_labelAlphabet.size(), m_options.rnnHiddenSize, m_options.hiddenSize);
  	m_classifier.setTagEmbFinetune(m_options.tagEmbFineTune);
  	m_classifier.setDropValue(m_options.dropProb);

  	//initial examples
	vector<vector<Example> > vecDevExamples(vecDevFiles.size()), vecTestExamples(vecTestFiles.size());
	for(int i = 0; i < vecDevInsts.size(); i++){
		initialExamples(vecDevInsts[i], vecDevExamples[i]);
		initialExamples(vecTestInsts[i], vecTestExamples[i]);
	}
	
	//set corpus weight size
	vector<int> vecSampleSize;
    str2int_vec(vecTrainCorpusSize, vecSampleSize);
	int inputSize = 0;
	for(int i = 0; i < vecSampleSize.size(); i++){
		if(vecSampleSize[i] < 1 || vecSampleSize[i] > vecTrainInstsPositions[i].size())
			vecSampleSize[i] = vecTrainInstsPositions[i].size();
		cout << "The " << i << "th Train Corpus Size:\t" << vecSampleSize[i] << endl; 
		inputSize += vecSampleSize[i];
	}
	
	cout << "All Train InputSize:\t" << inputSize << endl;
  	cout << endl;

	int batchBlock = inputSize / m_options.batchSize;
  	if (inputSize % m_options.batchSize != 0)
    		batchBlock++;

	//random shuffle
  	srand(0);
	vector<int> indexesAll;
	vector<vector<int> > indexesTrains(vecSampleSize.size());
  	
	for (int i = 0; i < inputSize; ++i)
    	indexesAll.push_back(i);
	for (int i = 0; i < vecSampleSize.size(); i++){
		for (int j = 0; j < vecTrainInstsPositions[i].size(); j++)
			indexesTrains[i].push_back(j);
	}

  	static Metric eval;
	static vector<Metric> metric_devs(vecDevFiles.size());
	static vector<Metric> metric_tests(vecTestFiles.size());//这边可能要修改
  	
	static vector<Example> allSampleTrainExamples, subExamples;
	static vector<Instance> allSampleTrainInstances;
	
	time_t start_iter_time, end_iter_time;

  	for (int iter = 0, stopIter = 0; iter < m_options.maxIter  && stopIter < m_options.stopIter; ++iter, ++stopIter) {
    	std::cout << "##### Iteration " << iter  << std::endl;
		time(&start_iter_time);
    	
		for(int i = 0; i < vecSampleSize.size(); i++){
			random_shuffle(indexesTrains[i].begin(), indexesTrains[i].end());	
		}
		random_shuffle(indexesAll.begin(), indexesAll.end());
        
      /*  cout << "After Shuffle:" << endl;
        for (int tmp_i =0; tmp_i < vecSampleSize.size(); tmp_i++){
            cout << "Train " << tmp_i << "th index: " ;
            for (int tmp_j =0; tmp_j < vecSampleSize[tmp_i]; tmp_j++){
                cout << indexesTrains[tmp_i][tmp_j] << "," ;
            }
            cout << endl;
        }
        cout << "Shuffle IndexesAll: ";
        for (int tmp_i=0; tmp_i < inputSize; tmp_i++){
            cout << indexesAll[tmp_i] << ",";
        }
        cout << endl;
        */
		
		allSampleTrainInstances.clear();
		allSampleTrainExamples.clear();
		
		size_t posi;
		for(int i_file = 0; i_file < vecTrainFiles.size(); i_file++){
			m_pipe.initInputFile(vecTrainFiles[i_file].c_str());
			for (int i_size = 0; i_size < vecSampleSize[i_file]; i_size++){
				posi = vecTrainInstsPositions[i_file][indexesTrains[i_file][i_size]];
				Instance *pInstance = m_pipe.readInstanceByPosi(posi);
				Instance tmp_trainInstance;
            	tmp_trainInstance.copyValuesFrom(*pInstance);
				allSampleTrainInstances.push_back(tmp_trainInstance);
			}
			m_pipe.uninitInputFile();
		}
		
		assert(allSampleTrainInstances.size() == inputSize);
		
		initialExamples(allSampleTrainInstances, allSampleTrainExamples);
		allSampleTrainInstances.clear();
    		
		eval.reset();

		//train update
    	for (int updateIter = 0; updateIter < batchBlock; updateIter++) {
      		subExamples.clear();
      		int start_pos = updateIter * m_options.batchSize;
      		int end_pos = (updateIter + 1) * m_options.batchSize;
      		if (end_pos > inputSize)
        		end_pos = inputSize;

      		for (int idy = start_pos; idy < end_pos; idy++) {
        		subExamples.push_back(allSampleTrainExamples[indexesAll[idy]]);
      		}

      		int curUpdateIter = iter * batchBlock + updateIter;
      		dtype cost = m_classifier.process(subExamples, curUpdateIter);

      		eval.overall_label_count += m_classifier._eval.overall_label_count;
      		eval.correct_label_count += m_classifier._eval.correct_label_count;

     		if ((curUpdateIter + 1) % m_options.verboseIter == 0) {
        		std::cout << "current: " << updateIter + 1 << ", total block: " << batchBlock << std::endl;
        		std::cout << "Cost = " << cost << ", Tag Correct(%) = " << eval.getAccuracy() << std::endl;
      		}	
			
			m_classifier.updateParams(m_options.belta1, m_options.belta2, m_options.regParameter, m_options.adaAlpha, m_options.adaEps);
    	}

		for (int i_dev = 0; i_dev < vecDevFiles.size(); i_dev++){
    		if (vecDevExamples[i_dev].size() > 0) {
      			bCurIterBetters[i_dev] = false;
      			if (!m_options.outBest.empty())
        			decodeInstResults.clear();
      			metric_devs[i_dev].reset();
      			for (int idx = 0; idx < vecDevExamples[i_dev].size(); idx++) {
        			vector<string> result_labels;
        			predict(vecDevExamples[i_dev][idx].m_features, result_labels, vecDevInsts[i_dev][idx].words);

        			if (m_options.seg)
          				vecDevInsts[i_dev][idx].SegEvaluate(result_labels, metric_devs[i_dev]);
        			else
          				vecDevInsts[i_dev][idx].Evaluate(result_labels, metric_devs[i_dev]);

        			if (!m_options.outBest.empty()) {
          				curDecodeInst.copyValuesFrom(vecDevInsts[i_dev][idx]);
          				curDecodeInst.assignLabel(result_labels);
          				decodeInstResults.push_back(curDecodeInst);
        			}
      			}	

				std::cout << vecDevFiles[i_dev] << " Dev:" << std::endl;
      			metric_devs[i_dev].print();

      			if (!m_options.outBest.empty() && metric_devs[i_dev].getAccuracy() > bestDISs[i_dev]) {
        			m_pipe.outputAllInstances(vecDevFiles[i_dev] + m_options.outBest, decodeInstResults);
        			bCurIterBetters[i_dev] = true;
      			}

      			if (vecTestExamples[i_dev].size() > 0) {
        			if (!m_options.outBest.empty())
          				decodeInstResults.clear();
        			metric_tests[i_dev].reset();
        			for (int idx = 0; idx < vecTestExamples[i_dev].size(); idx++) {
          				vector<string> result_labels;
          				predict(vecTestExamples[i_dev][idx].m_features, result_labels, vecTestInsts[i_dev][idx].words);

          				if (m_options.seg)
            				vecTestInsts[i_dev][idx].SegEvaluate(result_labels, metric_tests[i_dev]);
          				else
            				vecTestInsts[i_dev][idx].Evaluate(result_labels, metric_tests[i_dev]);

          				if (bCurIterBetters[i_dev] && !m_options.outBest.empty()) {
            				curDecodeInst.copyValuesFrom(vecTestInsts[i_dev][idx]);
            				curDecodeInst.assignLabel(result_labels);
            				decodeInstResults.push_back(curDecodeInst);
          				}
        			}
        			std::cout << vecTestFiles[i_dev] << " test:" << std::endl;
        			metric_tests[i_dev].print();

        			if (!m_options.outBest.empty() && bCurIterBetters[i_dev]) {
          				m_pipe.outputAllInstances(vecTestFiles[i_dev] + m_options.outBest, decodeInstResults);
        			}
      			}
      
      			if (m_options.saveIntermediate && metric_devs[i_dev].getAccuracy() > bestDISs[i_dev]) {
					stopIter = 0;
        			std::cout << vecDevFiles[i_dev] <<" Exceeds best previous performance of " << bestDISs[i_dev] << ". Saving news model file.." << std::endl;
        			bestDISs[i_dev] = metric_devs[i_dev].getAccuracy();
					stringstream ss;
					ss << iter;
					string s_updateIter = ss.str();
        			writeModelFile(modelFile + vecDevFiles[i_dev] +"_dev=" + s_updateIter);
      			}

    		}
		}

	time(&end_iter_time);
	cout << "one iter difftime:\t" << difftime(end_iter_time, start_iter_time) << "(s)" << endl;
  	}
}

int Labeler::predict(const vector<Feature>& features, vector<string>& outputs, const vector<string>& words) {
  assert(features.size() == words.size());
  vector<int> labelIdx, label2Idx;
  m_classifier.predict(features, labelIdx);
  outputs.clear();

  for (int idx = 0; idx < words.size(); idx++) {
    string label = m_labelAlphabet.from_id(labelIdx[idx]);
    outputs.push_back(label);
  }

  return 0;
}

void Labeler::test(const string& testFile, const string& outputFile, const string& modelFile) {
  loadModelFile(modelFile);
  vector<Instance> testInsts;
  m_pipe.readInstances(testFile, testInsts);

  vector<Example> testExamples;
  initialExamples(testInsts, testExamples);

  int testNum = testExamples.size();
  vector<Instance> testInstResults;
  Metric metric_test;
  metric_test.reset();
  for (int idx = 0; idx < testExamples.size(); idx++) {
    vector<string> result_labels;
    predict(testExamples[idx].m_features, result_labels, testInsts[idx].words);
    testInsts[idx].SegEvaluate(result_labels, metric_test);
    Instance curResultInst;
    curResultInst.copyValuesFrom(testInsts[idx]);
    testInstResults.push_back(curResultInst);
  }
  std::cout << "test:" << std::endl;
  metric_test.print();

  m_pipe.outputAllInstances(outputFile, testInstResults);

}

void Labeler::loadModelFile(const string& inputModelFile) {

}

void Labeler::writeModelFile(const string& outputModelFile) {

}

int main(int argc, char* argv[]) {
#if USE_CUDA==1
	InitTensorEngine();
#else
	InitTensorEngine<cpu>();
#endif

	std::string trainFiles = "", trainCorpusSize = ""; 
    std::string devFiles = "", testFiles = ""; 
	std::string modelFile = "";
	std::string optionFile = "";
	std::string outputFile = "";
	bool bTrain = false;
	dsr::Argument_helper ah;

	ah.new_flag("l", "learn", "train or test", bTrain);
	ah.new_named_string("train", "trainCorpus", "named_string",
			"training corpus to train a model, must when training", trainFiles);
	ah.new_named_string("trainSize", "trainCorpusSize", "named_string",
			"training corpus Size to train a model, must when training", trainCorpusSize);
	ah.new_named_string("dev", "devCorpus", "named_string",
			"development corpus to train a model, optional when training", devFiles);
	ah.new_named_string("test", "testCorpus", "named_string",
			"testing corpus to train a model or input file to test a model, optional when training and must when testing", testFiles);
	ah.new_named_string("model", "modelFile", "named_string",
			"model file, must when training and testing", modelFile);
	ah.new_named_string("option", "optionFile", "named_string",
			"option file to train a model, optional when training", optionFile);
	ah.new_named_string("output", "outputFile", "named_string",
			"output file to test, must when testing", outputFile);

	ah.process(argc, argv);

	Labeler tagger;
	if (bTrain) {
		tagger.train(trainFiles, trainCorpusSize, devFiles, testFiles, modelFile, optionFile);
	} else {
	}

	//test(argv);
	//ah.write_values(std::cout);
#if USE_CUDA==1
	ShutdownTensorEngine();
#else
	ShutdownTensorEngine<cpu>();
#endif
}
