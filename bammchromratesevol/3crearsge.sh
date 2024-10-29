#!/bin/bash
# My first script
mkdir sge
for bamm in {1..100} ; do
sed -e "s/bamm1/bamm$bamm/g" ./script_bamm.sge > ./sge/script_bamm$bamm.sge
done

