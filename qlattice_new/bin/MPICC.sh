#!/bin/bash
if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi
if [ intel = "$USE_COMPILER" ] ; then
    $run mpiicc "$@"
elif which mpiicc >/dev/null 2>&1 ; then
    if [ gcc = "$USE_COMPILER" ] ; then
        $run mpiicc -cc=gcc "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        $run mpiicc -cc=clang "$@"
    elif [ -f "/thfs1/home/fengxu/wxh/code/qlattice_new/bin/clang" ] ; then
        $run mpiicc -cc=clang "$@"
    elif [ -f "/thfs1/home/fengxu/wxh/code/qlattice_new/bin/gcc" ] ; then
        $run mpiicc -cc=gcc "$@"
    elif [ -f "/usr/bin/clang" ] ; then
        $run mpiicc -cc=clang "$@"
    elif [ -f "/usr/bin/gcc" ] ; then
        $run mpiicc -cc=gcc "$@"
    else
        $run mpiicc "$@"
    fi
else
    if [ gcc = "$USE_COMPILER" ] ; then
        OMPI_CC=gcc $run mpicc "$@"
    elif [ clang = "$USE_COMPILER" ] ; then
        OMPI_CC=clang $run mpicc "$@"
    else
        $run mpicc "$@"
    fi
fi
