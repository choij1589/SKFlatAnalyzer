#!/bin/bash
ERA=$1

SKFlat.py -a MeasTrigEff -i SingleElectron -n 10 -e 2016preVFP --userflags MeasMuLegs &
SKFlat.py -a MeasTrigEff -i SingleMuon -n 10 -e 2016preVFP --userflags MeasElLegs &
SKFlat.py -a MeasTrigEff -i MuonEG -n 10 -e 2016preVFP --userflags MeasEMuDZ &
SKFlat.py -a MeasTrigEff -i DoubleMuon -n 10 -e 2016preVFP --userflags MeasDblMuDZ &
SKFlat.py -a MeasTrigEff -l trigSamples.txt -n 20 -e 2016preVFP --userflags MeasMuLegs &
SKFlat.py -a MeasTrigEff -i trigSamples.txt -n 20 -e 2016preVFP --userflags MeasElLegs &
SKFlat.py -a MeasTrigEff -i trigSamples.txt -n 20 -e 2016preVFP --userflags MeasEMuDZ &
SKFlat.py -a MeasTrigEff -i trigSamples.txt -n 20 -e 2016preVFP --userflags MeasDblMuDZ &
