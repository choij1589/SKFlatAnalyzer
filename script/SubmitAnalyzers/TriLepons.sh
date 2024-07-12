#!/bin/bash
ERA=$1
CHANNEL=$2
MEMORY=$3

if [[ $CHANNEL == "Skim1E2Mu" ]]; then
    DATASTREAM="MuonEG"
elif [[ $CHANNEL == "Skim3Mu" ]]; then
    DATASTREAM="DoubleMuon"
else
    echo "Wrong channel $CHANNEL"
    exit 1
fi

SKFlat.py -a PromptEstimator -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
SKFlat.py -a MatrixEstimator -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst --memory $MEMORY --python &
SKFlat.py -a PromptEstimator -l SampleLists/triLepSamples.txt --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL,RunSyst, --memory $MEMORY --python &
