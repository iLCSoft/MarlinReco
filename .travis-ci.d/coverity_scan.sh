#!/bin/bash

ILCSOFT=/cvmfs/clicdp.cern.ch/iLCSoft/builds/current/CI_${COMPILER}
source $ILCSOFT/init_ilcsoft.sh

cd /Package
mkdir build
cd build
cmake -C $ILCSOFT/ILCSoft.cmake -DCMAKE_CXX_FLAGS="-fdiagnostics-color=always -DMARLINRECO_FORTRAN='OFF'" .. && \
export PATH=/Package/cov-analysis-linux64/bin:$PATH && \
cov-build --dir cov-int make -j2 && \
tar czvf myproject.tgz cov-int
