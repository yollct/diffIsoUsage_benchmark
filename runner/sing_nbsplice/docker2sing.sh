#!/bin/bash

docker build -t nbsplice ../docker_nbsplice/
docker run -v /var/run/docker.sock:/var/run/docker.sock -v /tmp/test:/output -t --privileged --rm quay.io/singularity/docker2singularity --name nbsplice nbsplice

cp /tmp/test/nbsplice.sif .