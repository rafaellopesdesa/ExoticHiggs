#!/bin/sh
if [ ! $?LD_LIBRARY_PATH ]; then
  export LD_LIBRARY_PATH=/Users/garren/lhc/LCG/hepmc-2.06.08/x86_64-mac106-gcc42-opt/lib
fi
if [ $?LD_LIBRARY_PATH ]; then
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/Users/garren/lhc/LCG/hepmc-2.06.08/x86_64-mac106-gcc42-opt/lib
fi
export PYTHIA8DATA=${PYTHIA8_HOME}/xmldoc
