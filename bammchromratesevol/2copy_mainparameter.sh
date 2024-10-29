#!/bin/bash
# My first script
for bamm in {1..100} ; do
sed -e "s/my_tree1.tree/my_tree$bamm.tree/g ; s/run_info1.txt/run_info$bamm.txt/g ; s/mcmc_out1.txt/mcmc_out$bamm.txt/g ; s/event_data1.txt/event_data$bamm.txt/g" ./Bamm.txt > ./"bamm$bamm"/Bamm.txt
done


#!/bin/bash
# My first script
for bamm in {1..100} ; do
sed -e "s/my_tree$bamm.tree/my_tree$bamm.tree/g"  ./my_trees/my_tree$bamm.tree > ./"bamm$bamm"/my_tree$bamm.tree
done

