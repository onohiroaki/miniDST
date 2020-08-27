#!/bin/sh

#
# setup environment to use right root, gcc, python and cmake for LCIO
#

. /cvmfs/sft.cern.ch/lcg/releases/gcc/8.2.0-3fa06/x86_64-slc6/setup.sh
export PATH=/cvmfs/ilc.desy.de/sw/x86_64_gcc82_sl6/CMake/3.15.5/bin/:$PATH
export PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-slc6-gcc8-opt/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_96/Python/2.7.16/x86_64-slc6-gcc8-opt/lib:${LD_LIBRARY_PATH}
. /cvmfs/ilc.desy.de/sw/x86_64_gcc82_sl6/root/6.18.04/bin/thisroot.sh
