#!/bin/bash
for cstart in $(seq -f%g 1 1 $1)
do
sbatch run.sh
done
