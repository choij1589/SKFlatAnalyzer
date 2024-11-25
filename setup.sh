#!/bin/bash
#export HOSTNAME=`hostname`
#if [[ $HOSTNAME == *"Mac"* ]]; then
#    RELEASE=""
#else
#HOSTNAME=`hostname`
RELEASE="`cat /etc/redhat-release`"
#fi

echo "@@@@ Working in $HOSTNAME"
if [[ $HOSTNAME == *"ai-tamsa"* ]]; then
  export SKFlat_WD="/data6/Users/$USER/SKFlatAnalyzer"
  export SKFlatRunlogDir="/gv0/Users/$USER/SKFlatRunlog"
  export SKFlatOutputDir="/gv0/Users/$USER/SKFlatOutput"
  # root configuration
  source ~/.conda-activate
  conda activate pyg
elif [[ $HOSTNAME == *"tamsa"* ]]; then
  export SKFlat_WD="/data6/Users/$USER/SKFlatAnalyzer"
  export SKFlatRunlogDir="/gv0/Users/$USER/SKFlatRunlog"
  export SKFlatOutputDir="/gv0/Users/$USER/SKFlatOutput"
  # root configuration
  # Singlarity image
  if [[ $RELEASE == *"Alma"* ]]; then
    source /opt/conda/bin/activate
    conda activate pyg
  else
    # temporarily use ROOT and python from LCG environment
    source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc12-opt/setup.sh
  fi
elif [[ $HOSTNAME == *"cms"* ]]; then
  export SKFlat_WD="/data6/Users/$USER/SKFlatAnalyzer"
  export SKFlatRunlogDir="/data6/Users/$USER/SKFlatRunlog"
  export SKFlatOutputDir="/data6/Users/$USER/SKFlatOutput"
  # root configuration
  source /home/choij/miniconda3/bin/activate
  conda activate pyg
elif [[ $HOSTNAME == *"private"* ]]; then
  export SKFlat_WD="/home/$USER/workspace/SKFlatAnalyzer"
  export SKFlatRunlogDir="/home/$USER/workspace/SKFlatRunlog"
  export SKFlatOutputDir="/home/$USER/workspace/SKFlatOutput"
  # root configuration
  source ~/.conda-activate
  conda activate pyg
elif [[ $HOSTNAME == *"Mac"* ]]; then
  export SKFlat_WD="/Users/$USER/workspace/SKFlatAnalyzer"
  export SKFlatRunlogDir="/Users/$USER/workspace/SKFlatRunlog"
  export SKFlatOutputDir="/Users/$USER/workspace/SKFlatOutput"
  # root configuration
  source ~/myenv/bin/activate
  export ROOTSYS=$(brew --prefix root)
  export PATH=$ROOTSYS/bin:$PATH
  export PYTHONPATH=$ROOTSYS/lib/root:$PYTHONPATH
else
  echo "Unknown hostname $HOSTNAME"
  exit 1
fi

export SKFlat_LIB_PATH=$SKFlat_WD/lib/
mkdir -p $SKFlat_LIB_PATH
mkdir -p $SKFlat_WD/tar

export SKFlatV="Run2UltraLegacy_v3"
mkdir -p $SKFlat_WD/data/$SKFlatV
export DATA_DIR=$SKFlat_WD/data/$SKFlatV

alias skout="cd $SKFlatOutputDir/$SKFlatV"
export MYBIN=$SKFlat_WD/bin/
export PYTHONDIR=$SKFlat_WD/python/
export PYTHONPATH="${PYTHONPATH}:${PYTHONDIR}"

# setting LHAPDF
if [[ -d "external/lhapdf" ]]; then
    #export PATH=$PATH:$SKFlat_WD/external/lhapdf/bin
    #export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKFlat_WD/external/lhapdf/lib
    #export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKFlat_WD/external/lhapdf/lib
    export LHAPDFDIR=$SKFlat_WD/external/lhapdf
    export LHAPDF_DATA_PATH=$SKFlat_WD/external/lhapdf/share/LHAPDF
else
    export LHAPDFDIR=$SKFlat_WD/external/lhapdf
    export LHAPDF_DATA_PATH=$LHAPDFDIR/data
fi
export LHAPDF_INCLUDE_DIR=`lhapdf-config --incdir`
export LHAPDF_LIB_DIR=`lhapdf-config --libdir`

echo "@@@@ LHAPDF include: $LHAPDF_INCLUDE_DIR"
echo "@@@@ LHAPDF lib: $LHAPDF_LIB_DIR"
echo "@@@@ reading data from $LHAPDF_DATA_PATH"

export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$SKFlat_WD/DataFormats/include/:$SKFlat_WD/AnalyzerTools/include/:$SKFlat_WD/Analyzers/include/
export PATH=${MYBIN}:${PYTHONDIR}:${LHAPDFDIR}/bin:${PATH}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SKFlat_LIB_PATH:$LHAPDF_LIB_DIR
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$SKFlat_LIB_PATH:$LHAPDF_LIB_DIR
source $SKFlat_WD/bin/BashColorSets.sh
## submodules ##
#source bin/CheckSubmodules.sh

if [ "$1" = "-q" ];then
    return
fi

## Todo list ##
python python/PrintToDoLists.py
source $SKFlat_WD/tmp/ToDoLists.sh
rm $SKFlat_WD/tmp/ToDoLists.sh

CurrentGitBranch=`git branch | grep \* | cut -d ' ' -f2`
printf "> Current SKFlatAnalyzer branch : "${BRed}$CurrentGitBranch${Color_Off}"\n"
echo "-----------------------------------------------------------------"
## Log Dir ##
# echo "* Your Log Directory Usage (ctrl+c to skip)"
# du -sh $SKFlatRunlogDir
