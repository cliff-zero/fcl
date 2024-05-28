#!/bin/bash
if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi
if [ intel = "$USE_COMPILER" ] ; then
    $run mpiicpc "$@"
elif which mpiicpc >/dev/null 2>&1 ; then
    if [ gcc = "$USE_COMPILER" ] ; then
        $run mpiicpc -cxx=g++ "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        $run mpiicpc -cxx=clang++ "$@"
    elif [ -f "/thfs1/home/fengxu/wxh/code/qlattice_new/bin/clang++" ] ; then
        $run mpiicpc -cxx=clang++ "$@"
    elif [ -f "/thfs1/home/fengxu/wxh/code/qlattice_new/bin/g++" ] ; then
        $run mpiicpc -cxx=g++ "$@"
    elif [ -f "/usr/bin/clang++" ] ; then
        $run mpiicpc -cxx=clang++ "$@"
    elif [ -f "/usr/bin/g++" ] ; then
        $run mpiicpc -cxx=g++ "$@"
    else
        $run mpiicpc "$@"
    fi
elif which mpicxx >/dev/null 2>&1 ; then
    if [ gcc = "$USE_COMPILER" ] ; then
        OMPI_CXX=g++ $run mpicxx "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        OMPI_CXX=clang++ $run mpicxx "$@"
    else
        $run mpicxx "$@"
    fi
else
    if [ gcc = "$USE_COMPILER" ] ; then
        OMPI_CXX=g++ $run mpic++ "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        OMPI_CXX=clang++ $run mpic++ "$@"
    else
        $run mpic++ "$@"
    fi
fi
