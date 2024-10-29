#!/bin/bash
# My first script
#mkdir sge
for bamm in {1..100} ; do
sed -e "s/bamm1/bamm$bamm/g; s/Script1.R/Script$bamm.R/g" ./script_Rbamm.sge > ./sge/script_Rbamm$bamm.sge
done

