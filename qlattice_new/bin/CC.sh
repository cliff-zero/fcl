#!/bin/bash
if which python3 >/dev/null 2>&1 ; then
    run=compiler-options.py
else
    run=
fi
if [ intel = "$USE_COMPILER" ] ; then
    $run icc "$@"
elif [ gcc = "$USE_COMPILER" ] ; then
    $run gcc "$@"
elif [ clang = "$USE_COMPILER" ] ; then
    $run clang "$@"
elif [ -f "/thfs1/home/fengxu/wxh/code/qlattice_new/bin/clang" ] ; then
    $run clang "$@"
elif [ -f "/thfs1/home/fengxu/wxh/code/qlattice_new/bin/gcc" ] ; then
    $run gcc "$@"
elif which icc >/dev/null 2>&1 ; then
    $run icc "$@"
elif [ -f "/usr/bin/clang" ] ; then
    $run clang "$@"
else
    $run gcc "$@"
fi
