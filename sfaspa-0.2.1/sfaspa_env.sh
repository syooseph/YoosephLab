#!/bin/bash

ulimit -s unlimited

if [[ $SPA_HOME == "" ]]; then 
	echo "SPA_HOME is not properly configured"
	exit
fi

export SPA_BIN=$SPA_HOME/bin
export SPA_SCRIPT=$SPA_HOME/script
export FGS=$SPA_HOME/3rd/FragGeneScan1.16
export MGA=$SPA_HOME/3rd/mga_ia64
export PATH=$SPA_BIN:$SPA_SCRIPT:$FGS:$MGA:$PATH

export LD_LIBRARY_PATH=$SPA_HOME/3rd/libdivsufsort-2.0.1/build/lib:$LD_LIBRARY_PATH
