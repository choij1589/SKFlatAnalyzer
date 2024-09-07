#!/bin/bash
ERA=$1
CHANNEL=$2

if [[ $CHANNEL == "RunEMu" ]]; then
    DATASTREAM="MuonEG"
elif [[ $CHANNEL == "RunDiMu" ]]; then
    DATASTREAM="DoubleMuon"
else
    echo "Wrong channel $CHANNEL"
    exit 1
fi

SKFlat.py -a CR_DiLepton -i $DATASTREAM -n 20 -e ${ERA} --userflags $CHANNEL --python &
SKFlat.py -a CR_DiLepton -l SampleLists/diLepSamples.txt -n 20 -e ${ERA} --userflags $CHANNEL,RunSyst --python &
