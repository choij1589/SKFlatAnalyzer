#!/bin/bash
ERA=$1
CHANNEL=$2
DATASTREAM=""
QCD=""

if [[ $ERA == "2016"* ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
        QCD="qcdMuEnriched"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="DoubleEG"
        QCD="qcdEMEnriched"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV2 -i $DATASTREAM -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i TTLJ_powheg -n 30 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i TTLL_powheg -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WJets_MG    -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i DYJets      -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i DYJets10to50_MG   -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WW_pythia         -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WZ_pythia         -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i ZZ_pythia         -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_sch_Lep -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tch_top_Incl          -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tch_antitop_Incl      -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tW_top_NoFullyHad     -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tW_antitop_NoFullyHad -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -l SampleLists/${QCD}.txt          -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
elif [[ $ERA == "2017" ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
        QCD="qcdMuEnriched"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="SingleElectron"
        QCD="qcdEMEnriched"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV2 -i $DATASTREAM -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i TTLJ_powheg -n 30 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i TTLL_powheg -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WJets_MG    -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i DYJets      -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i DYJets10to50_MG   -n 10 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WW_pythia         -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WZ_pythia         -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i ZZ_pythia         -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_sch_Lep -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tch_top_Incl          -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tch_antitop_Incl      -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tW_top_NoFullyHad     -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tW_antitop_NoFullyHad -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -l SampleLists/${QCD}.txt          -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
elif [[ $ERA == "2018" ]]; then
    if [[ $CHANNEL == "MeasFakeMu"* ]]; then
        DATASTREAM="DoubleMuon"
        QCD="qcdMuEnriched"
    elif [[ $CHANNEL == "MeasFakeEl"* ]]; then
        DATASTREAM="EGamma"
        QCD="qcdEMEnriched"
    else
        echo "Wrong channel $CHANNEL"
        exit 1
    fi
    SKFlat.py -a MeasFakeRateV2 -i $DATASTREAM -n 15 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i TTLJ_powheg -n 40 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i TTLL_powheg -n 30 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WJets_MG    -n 30 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i DYJets      -n 30 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i DYJets10to50_MG   -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WW_pythia         -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i WZ_pythia         -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i ZZ_pythia         -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_sch_Lep -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tch_top_Incl          -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tch_antitop_Incl      -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tW_top_NoFullyHad     -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -i SingleTop_tW_antitop_NoFullyHad -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
    SKFlat.py -a MeasFakeRateV2 -l SampleLists/${QCD}.txt          -n 20 -e $ERA --userflags $CHANNEL,RunSystSimple &
else
    echo "Wrong era $ERA"
    exit 1
fi
