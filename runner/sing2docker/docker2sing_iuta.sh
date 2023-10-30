#!/bin/bash

#docker build -t iuta ../docker_nbsplice/
docker run -v /var/run/docker.sock:/var/run/docker.sock -v /tmp/test:/output -t --privileged --rm quay.io/singularity/docker2singularity --name iuta iuta

cp /tmp/test/iuta.sif .