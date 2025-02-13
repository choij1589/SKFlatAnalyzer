# SKFlatAnalyzer

## Manual

https://jskim.web.cern.ch/jskim/SKFlat/Manual/Manual_SKFlat.pdf

## TODO
- To reduce CPU overload, update SKFlat.py to use condor\_wait rather than use python to read the condor outputs

## Where to put the analyzer
TAMSA 1/2 : /data6/Users/$USER/
Only supports for TAMSA server now.

## First time setup
Before clone the repository, it is a good strategy to fork from the upstream repository.
```bash
#### When first time git clone, use the option "--recursive" to initiate the submodules
git clone --recursive git@github.com:<gitaccount>/SKFlatAnalyzer.git
cd SKFlatAnalyzer

#### First time setup script
source bin/FirstTimeSetup.sh 
source setup.sh

#### Install lhapdf if needed
./script/build_lhapdf.sh

#### You have to edit user info
#### First, copy the temply using the command below
cp $SKFlat_WD/python/UserInfo_template.py $SKFlat_WD/python/UserInfo_${USER}.py 
#### Then, edit $SKFlat_WD/python/UserInfo_${USER}.py
```
Compile
> Note that after submitting condor jobs, we use singularity image based on conda setup.
> Compiliation should be done inside the singularity image to match with the environments.

```bash
#### Basic build
./script/build.sh

#### Using singularity image
singularity exec /data9/Users/choij/Singularity/images/cuda11.8 script/build.sh

#### Now, run setup script.
#### This should be done for every new shell
source setup.sh
```

## Test job
```bash
# C++ based analyzers
SKFlat.py -a ExampleRun -i DYJets -n 50 -e 2017 &
# python based analyzers
SKFlat.py -a TutorialRun -i DYJets -n 50 -e 2017 --python &
```

## Making a new Ananlyzer
```bash
cd script/MakeCycleSkeleton/
```
Then, run
```bash
python MakeCycleSkeleton.py NewAnalyzer # <- put new analyzer name
```
It will print below lines (execute the lines) :
```bash
mv NewAnalyzer.h $SKFlat_WD/Analyzers/include/
mv NewAnalyzer.C $SKFlat_WD/Analyzers/src/
```

Then, add
```bash
#pragma link C++ class NewAnalyzer+;
```
in Analyzers/include/Analyzers_LinkDef.h

For python based analyzers, check
```bash
PyAnalyzers/TutorialRun.py
```

## Detailed descriptions

Look Analyzers/src/ExampleRun.C

## Adding samples
To add a sample, you should add two files.  
```bash
data/$SKFlatV/$ERA/Sample/ForSNU/$ALIAS.txt               # list of file paths.  
data/$SKFlatV/$ERA/Sample/CommonSampleInfo/$ALIAS.txt     # alias, DAS name, cross section, nevent, sum(sign) and sum(weight).
```
And one file should be edited
```bash
data/$SKFlatV/$ERA/Sample/SampleSummary_*.txt:  # CommonSampleInfo files in one file. This file is not actually used by SKFlatAnalyzer. It is just a summary for users.  
```
You can do it manually, or use scripts as below. The scripts use SampleSummary files for all SKFlat versions to find alias and cross section. If the sample is never used before the scripts will ask you alias and cross section. If it is annoying you can make temporal SampleSummary file like ```data/$SKFlatV/$ERA/Sample/SampleSummary_temp.txt``` which containing alias, DAS name and cross section (other information is not needed).

1. Make the file list file using bin/UpdateSampleForSNU.sh script.
```bash
./bin/UpdateSampleForSNU.sh $SAMPLEDIRECTORY
```
2. Run GetEffLumi analyzer to get nevent, sum(sign) and sum(weight)
```bash
SKFlat.py -a GetEffLumi -e $ERA -n 20 -i $SAMPLEALIAS
```
3. Make CommonSampleInfo file using bin/UpdateCommonSampleInfo.sh script. It reads the root file outputs ran with singularity image, so it would be good to used singularity to avoid any ROOT version mismatch.
```bash
# Renew your shell!!!
singularity shell --bind /gv0/Users/choij/SKFlatOutput /data9/Users/choij/Singularity/images/cpuonly
source setup.sh
./bin/UpdateCommonSampleInfo.sh
```
4. Update SampleSummary file using Summarize.py
```bash
cd data/$SKFlatV/$ERA/Sample
python Summarize.py
```
# Tips
## Skimming samples
It is a good idea to first skim the sample in dedicated directory and move it /gv0/DATA/SKFlat, e.g.:
```bash
SKFlat.py -a SkimTree_SS2lOR3l -i DYJets -e 2017 -n 30 -o /gv0/Users/choij/temp
```

## Making PR
Start from the CMSSNU's master branch of CMSSNU when making pull request to prevent anoying conflicts.

