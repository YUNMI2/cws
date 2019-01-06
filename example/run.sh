#!/bin/bash 
#absolute path point to NNSegmentation directory
workspace=`pwd`
tooldir=/data/yzhu/Experiment/WS-Domain-Adaptation/Add-News/ZhuXian/CTB5+BK+ZX-sparse/
corpus=ZhuXian
outputdir=$corpus.sample

mkdir -p $outputdir
rm $outputdir/* -rf

function extract
{
    #extracting your features here
    echo "[self implementation]"
}

function runLSTM
{
    cmd=$1
    option=$2

    mkdir $workspace/$outputdir/$cmd.$corpus -p
    ln -s $workspace/sparsefeats/CTB50T.sparse.bichar.feats $workspace/$outputdir/$cmd.$corpus/ctb5.train.feats 
    ln -s $workspace/sparsefeats/baike.zxd.notrainwords.nomoretrainchars.nobou.sparse.bichar.feats $workspace/$outputdir/$cmd.$corpus/bk.train.feats
    ln -s $workspace/sparsefeats/zhuxian.ctb.raw.mm.sentences $workspace/$outputdir/$cmd.$corpus/zx.train.feats 
    ln -s $workspace/sparsefeats/ZXD.sparse.bichar.feats $workspace/$outputdir/$cmd.$corpus/ZXD.sparse.bichar.feats
    ln -s $workspace/sparsefeats/ZXE.sparse.bichar.feats $workspace/$outputdir/$cmd.$corpus/ZXE.sparse.bichar.feats
    echo $tooldir/$cmd 
	echo $workspace/$outputdir/$cmd/
	cp $tooldir/$cmd $workspace/$outputdir/$cmd.$corpus/
    	train_file=$workspace/$outputdir/$cmd.$corpus/ctb5.train.feats,$workspace/$outputdir/$cmd.$corpus/bk.train.feats,$workspace/$outputdir/$cmd.$corpus/zx.train.feats
        train_corpus=10000,10000,10000
    	dev_file=$workspace/$outputdir/$cmd.$corpus/ZXD.sparse.bichar.feats
    	test_file=$workspace/$outputdir/$cmd.$corpus/ZXE.sparse.bichar.feats
    #character bigram embedding and character trigram embedding should use a comma to separate 
    nohup $workspace/$outputdir/$cmd.$corpus/$cmd -l \
        -train $train_file \
        -trainSize $train_corpus \
        -dev $dev_file \
        -test $test_file \
        -option $option \
        -model $workspace/$outputdir/$cmd/$cmd.model \
    >$workspace/$outputdir/$cmd.log 2>&1 &
}

echo "Step 1: Extracting Features..."
extract $corpus 
echo "Step 2: Running LSTMCRFMLLabeler..."
cmds="LSTMCRFMLLabeler"
runLSTM LSTMCRFMLLabeler ./options/option.lstm
#runLSTM SparseLSTMCRFMMLabeler ./options/option.sparse+lstm
echo "Successfully run!"
