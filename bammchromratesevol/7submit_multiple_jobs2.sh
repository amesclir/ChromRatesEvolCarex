#!/bin/bash
for bamm in {1..100}; do
qsub bamm/bammchromratesevol/sge/script_Rbamm$bamm.sge; done
