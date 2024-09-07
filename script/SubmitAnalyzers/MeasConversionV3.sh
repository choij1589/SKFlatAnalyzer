#!/bin/bash
ERA=$1
CHANNEL=$2

if [[ $CHANNEL == "Skim1E2Mu" ]]; then
    DATASTREAM="MuonEG"
elif [[ $CHANNEL == "Skim3Mu" ]]; then
    DATASTREAM="DoubleMuon"
else
    echo "Wrong channel $CHANNEL"
    exit 1
fi
SKFlat.py -a MeasConversionV3 -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConvMatrixV3 -i $DATASTREAM --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL  --python &
SKFlat.py -a MeasConversionV3 -l SampleLists/triLepSamples.txt --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags $CHANNEL --python &
