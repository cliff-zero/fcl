#!/bin/bash
if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi
if [ intel = "$USE_COMPILER" ] ; then
    $run icpc "$@"
elif [ gcc = "$USE_COMPILER" ] ; then
    $run g++ "$@"
elif [ clang = "$USE_COMPILER" ] ; then
    $run clang++ "$@"
elif [ -f "/thfs1/home/fengxu/wxh/code/qlattice_new/bin/clang++" ] ; then
    $run clang++ "$@"
elif [ -f "/thfs1/home/fengxu/wxh/code/qlattice_new/bin/g++" ] ; then
    $run g++ "$@"
elif which icpc >/dev/null 2>&1 ; then
    $run icpc "$@"
elif [ -f "/usr/bin/clang++" ] ; then
    $run clang++ "$@"
else
    $run g++ "$@"
fi
