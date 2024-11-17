#!/bin/bash
ERA=$1
CHANNEL=$2
MEMORY=2000

if [[ $CHANNEL == "Skim1E2Mu" ]]; then
    DATASTREAM="MuonEG"
elif [[ $CHANNEL == "Skim3Mu" ]]; then
    DATASTREAM="DoubleMuon"
else
    echo "Wrong channel $CHANNEL"
    exit 1
fi

SKFlat.py -a PromptSkimmer -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
SKFlat.py -a MatrixSkimmer -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
SKFlat.py -a PromptSkimmer -l SampleLists/triLepSamples.txt --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --memory $MEMORY --python &
SKFlat.py -a PromptSkimmer -l SampleLists/signalSamples.txt -n 10 -e ${ERA} --userflags $CHANNEL,RunTheoryUnc --memory $MEMORY --python &
