FROM python:2
FROM ubuntu:18.04

# RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools unar libbz2-dev liblzma-dev libglpk-dev libnetcdf-c++4-1 git libfontconfig1-dev cmake octave-dev 
RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y wget zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev libncurses5-dev samtools libssl1.0.0 cmake make
RUN apt-get --no-install-recommends --fix-broken install -y octave liboctave-dev
RUN apt-get --no-install-recommends --fix-broken install -y cruft

ENV CONDA_DIR /opt/conda
RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \ 
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda install -c bioconda samtools
RUN conda install -c conda-forge libdeflate

RUN ln -s /opt/conda/pkgs/libdeflate-1.17-h5eee18b_0/bin/libdeflate.so.0 /lib/x86_64-linux-gnu/libdeflate.so.0
COPY ./rDiff /rDiff
#RUN apt update && apt install parallel -y --force-yes
WORKDIR /rDiff
RUN pip install --no-cache-dir numpy scipy 
# Installing Mamba
COPY . . 

WORKDIR /MOUNT
#ENTRYPOINT SCRIPT
COPY --chown=nobody:nogroup ./ENTRYPOINT.sh /ENTRYPOINT.sh


RUN chmod 777 /ENTRYPOINT.sh
ENTRYPOINT ["/ENTRYPOINT.sh"]
