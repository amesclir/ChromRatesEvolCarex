#!/bin/bash
# My first script
for bamm in {1..100} ; do
sed -e "s/my_tree1.tree/my_tree$bamm.tree/g ; s/mcmc_out1.txt/mcmc_out$bamm.txt/g ; s/event_data1.txt/event_data$bamm.txt/g ; s/branch_matrix1.csv/branch_matrix$bamm.csv/g" ./Script.R > ./bamm$bamm/Script$bamm.R
done

