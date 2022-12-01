#!/bin/bash

docker run -v /var/run/docker.sock:/var/run/docker.sock -v /tmp/test:/output -t --privileged --rm quay.io/singularity/docker2singularity --name jcseq jcseq

cp /tmp/test/jcseq.sif .