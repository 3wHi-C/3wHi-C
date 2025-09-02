#!/bin/bash

for i in {1..22} X
do
	/usr/local/bin/Rscript 3w_call_TAD_boundary_IC1.R chr${i} &
	/usr/local/bin/Rscript 3w_calculate_boundary_delta_0.1_IC1.R chr${i} &
done
