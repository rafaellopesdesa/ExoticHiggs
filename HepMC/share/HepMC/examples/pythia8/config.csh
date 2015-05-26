#!/bin/csh
if( ! $?LD_LIBRARY_PATH ) then
  setenv LD_LIBRARY_PATH /Users/garren/lhc/LCG/hepmc-2.06.08/x86_64-mac106-gcc42-opt/lib
else
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/Users/garren/lhc/LCG/hepmc-2.06.08/x86_64-mac106-gcc42-opt/lib
endif
setenv PYTHIA8DATA ${PYTHIA8_HOME}/xmldoc
