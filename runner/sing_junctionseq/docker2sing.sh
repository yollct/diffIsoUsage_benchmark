#!/bin/bash

docker run -v /var/run/docker.sock:/var/run/docker.sock -v /tmp/test:/output -w --privileged --rm quay.io/singularity/docker2singularity -m "/nfs/scratch/chit/isbench_covid/" --name jcseq jcseq

cp /tmp/test/jcseq.sif .

