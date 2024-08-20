#!/bin/bash
ERA=$1

if [ $ERA = "2018" ]; then
    SKFlat.py -a MeasTrigEff -i EGamma -n 10 -e $ERA --userflags MeasMuLegs &
else
    SKFlat.py -a MeasTrigEff -i SingleElectron -n 10 -e $ERA --userflags MeasMuLegs &
fi
SKFlat.py -a MeasTrigEff -i SingleMuon -n 10 -e $ERA --userflags MeasElLegs &
SKFlat.py -a MeasTrigEff -i MuonEG -n 10 -e $ERA --userflags MeasEMuDZ &
SKFlat.py -a MeasTrigEff -i DoubleMuon -n 10 -e $ERA --userflags MeasDblMuDZ &
SKFlat.py -a MeasTrigEff -l SampleLists/trigSamples.txt -n 20 -e $ERA --userflags MeasMuLegs &
SKFlat.py -a MeasTrigEff -l SampleLists/trigSamples.txt -n 20 -e $ERA --userflags MeasElLegs &
SKFlat.py -a MeasTrigEff -l SampleLists/trigSamples.txt -n 20 -e $ERA --userflags MeasEMuDZ &
SKFlat.py -a MeasTrigEff -l SampleLists/trigSamples.txt -n 20 -e $ERA --userflags MeasDblMuDZ &
