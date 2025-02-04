#!/bin/bash
ERA=$1

SKFlat.py -a JetIDandTaggingStudy -l SampleLists/signalSamples.txt -n 10 -e ${ERA} --userflags Skim1E2Mu --python &
SKFlat.py -a JetIDandTaggingStudy -i TTLL_powheg --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags Skim1E2Mu --python &
SKFlat.py -a JetIDandTaggingStudy -i DYJets --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags Skim1E2Mu --python &
SKFlat.py -a JetIDandTaggingStudy -i ttZToLLNuNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags Skim1E2Mu --python &

SKFlat.py -a JetIDandTaggingStudy -l SampleLists/signalSamples.txt -n 10 -e ${ERA} --userflags Skim3Mu --python &
SKFlat.py -a JetIDandTaggingStudy -i TTLL_powheg --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags Skim3Mu --python &
SKFlat.py -a JetIDandTaggingStudy -i DYJets --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags Skim3Mu --python &
SKFlat.py -a JetIDandTaggingStudy -i ttZToLLNuNu --skim SkimTree_SS2lOR3l -n 10 -e ${ERA} --userflags Skim3Mu --python &
