#!/bin/bash
ERAs=("2016preVFP" "2016postVFP" "2017" "2018")
CHANNELs=("Skim1E2Mu" "Skim3Mu")

for ERA in ${ERAs[@]}; do
    SKFlat.py -a AcceptanceStudy -l SampleLists/signalSamples.txt -e $ERA -n 10 --userflags Skim1E2Mu &
    SKFlat.py -a AcceptanceStudy -l SampleLists/signalSamples.txt -e $ERA -n 10 --userflags Skim3Mu &
done
