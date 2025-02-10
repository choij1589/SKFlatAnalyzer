#!/bin/bash
# clear the lhapdf directory
rm -rf $SKFlat_WD/external/lhapdf

# download and install lhapdf
cd $SKFlat_WD/external
wget https://lhapdf.hepforge.org/downloads/\?f\=LHAPDF-6.5.5.tar.gz -O LHAPDF-6.5.5.tar.gz
tar xf LHAPDF-6.5.5.tar.gz
cd LHAPDF-6.5.5
ARCH=`arch`
if [ $ARCH == "arm64" ]; then
    ./configure --prefix=$SKFlat_WD/external/lhapdf/osx CXX=clang++
else
    ./configure --prefix=$SKFlat_WD/external/lhapdf/redhat CXX=g++
fi
make -j8 && make install
cd $SKFlat_WD/external
rm -rf LHAPDF*

# download CMS central PDF sets
mkdir -p $SKFlat_WD/external/lhapdf/data
cd $SKFlat_WD/external/lhapdf/data
wget https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_hessian_pdfas.tar.gz
wget https://lhapdfsets.web.cern.ch/current/NNPDF31_nlo_hessian_pdfas.tar.gz
tar xf NNPDF31_nnlo_hessian_pdfas.tar.gz
tar xf NNPDF31_nlo_hessian_pdfas.tar.gz
cd $SKFlat_WD
