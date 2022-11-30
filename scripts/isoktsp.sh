#!/bin/bash



java -jar ../../iso-kTSP_v1.0.3.jar ../simulated_reads/isoktsp_count.txt -i -n 2 -s 5000 -o ../results/salmon_isoktsp_output.txt --seed 1234

java -jar ../../iso-kTSP_v1.0.3.jar ../simulated_reads/isoktsp_count.txt -i -n 2 -s 5000 -o ../results/rsem_isoktsp_output.txt --seed 1234