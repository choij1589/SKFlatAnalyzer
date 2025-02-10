#!/bin/bash
# check os
if [[ "$(uname)" == "Darwin" ]]; then
    export SYSTEM="osx"
elif [[ -f "/etc/redhat-release" ]]; then
    export SYSTEM="redhat"
else
    echo "Unsupported OS"
    return 1
fi

HOSTNAME=`hostname`
RELEASE="`cat /etc/redhat-release`"

echo "@@@@ Working in $HOSTNAME"
if [[ $HOSTNAME == *"tamsa"* ]]; then
  export SKFlat_WD="/data9/Users/$USER/Sync/workspace/SKFlatAnalyzer"
  export SKFlatRunlogDir="/gv0/Users/$USER/SKFlatRunlog"
  export SKFlatOutputDir="/gv0/Users/$USER/SKFlatOutput"
  # root configuration
  if [[ -n "$APPTAINER_NAME" || -n "$SINGULARITY_NAME" ]]; then
    export PATH="/opt/conda/bin:${PATH}"
    export MAMBA_ROOT_PREFIX="/opt/conda"
    eval "$(micromamba shell hook -s zsh)"
    micromamba activate Nano
  else
    export PATH="$HOME/micromamba/bin:${PATH}"
    export MAMBA_ROOT_PREFIX="$HOME/micromamba"
    eval "$(micromamba shell hook -s zsh)"
    micromamba activate Nano
    # temporarily use ROOT and python from LCG environment
    #source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc12-opt/setup.sh
  fi
elif [[ $HOSTNAME == *"private"* ]]; then
  export SKFlat_WD="/home/$USER/Sync/workspace/SKFlatAnalyzer"
  export SKFlatRunlogDir="/home/$USER/workspace/SKFlatRunlog"
  export SKFlatOutputDir="/home/$USER/workspace/SKFlatOutput"
  # root configuration
  export PATH="$HOME/micromamba/bin:${PATH}"
  export MAMBA_ROOT_PREFIX="$HOME/micromamba"
  eval "$(micromamba shell hook -s zsh)"
  micromamba activate Nano
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
export PATH=$PATH:$SKFlat_WD/external/lhapdf/$SYSTEM/bin
export LHAPDFDIR=$SKFlat_WD/external/lhapdf/$SYSTEM
export LHAPDF_DATA_PATH=$SKFlat_WD/external/lhapdf/data
export LHAPDF_INCLUDE_DIR=$SKFlat_WD/external/lhapdf/$SYSTEM/include
export LHAPDF_LIB_DIR=$SKFlat_WD/external/lhapdf/$SYSTEM/lib

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
