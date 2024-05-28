#!/bin/bash
for jid in $(seq $1 1 $2)
do
yhcancel ${jid}
done
