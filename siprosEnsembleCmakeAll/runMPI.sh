#!/bin/bash

currentPath=$(pwd)
echo $currentPath
# module load CMake/3.16.4-GCCcore-9.3.0 OpenMPI/4.0.3-GCC-9.3.0
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${currentPath}/libSiprosEnsembleMPI
echo $LD_LIBRARY_PATH
mpirun ./SiprosEnsembleMPI "$@"
