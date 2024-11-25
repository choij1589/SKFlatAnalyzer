#!/bin/bash
if [[ "$HOSTNAME" == "tamsa"* ]]; then
    source /data6/Users/choij/SKFlatAnalyzer/setup.sh
fi

export SKFlat_BUILDDIR=$SKFlat_WD/build

echo @@@@ Prepare to build SKFlatAnalyzer in $SKFlat_BUILDDIR
rm -rf $SKFlat_BUILDDIR $SKFlat_LIB_PATH
mkdir $SKFlat_BUILDDIR && cd $SKFlat_BUILDDIR

cmake -DCMAKE_INSTALL_PREFIX=$SKFlat_WD $SKFlat_WD

echo @@@@ make -j6
make -j6

echo @@@@ install to $SKFlat_LIB_PATH
make install
