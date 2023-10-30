FROM rocker/r-ver:3.3.3

# FROM condaforge/mambaforge:4.9.2-5 as conda

# RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y libcurl4-openssl-dev libssl-dev libxml2-dev 
# RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y make curl

COPY ./JunctionSeq/JunctionSeq /JunctionSeq

# RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \ 
#     && bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda \
#     && rm -f Miniconda3-latest-Linux-x86_64.sh 
# RUN echo "R_LIBS_USER='/Rlib'" >> /usr/local/lib/R/etc/Renviron
# RUN echo "R_LIBS='/Rlib'" >> /usr/local/lib/R/etc/Renviron

# RUN R -e "install.packages('plyr', version='1.8.3', dependencies=TRUE, repos='https://cloud.r-project.org')"  
# RUN R -e "install.packages('doParallel', version='1.0.10', dependencies=TRUE, repos='https://cloud.r-project.org')"  
# RUN R -e "install.packages('foreach', version='1.4.3', dependencies=TRUE, repos='https://cloud.r-project.org')"
# RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('GenomicRanges')" 
# RUN R -e "source('http://bioconductor.org/biocLite.R'); biocLite('DESeq2')" 

# RUN conda create -n r3.3 -c conda-forge r-base=3.3.2 && \
#     echo "conda activate r3.3" >> ~/.bashrc 

#RUN conda install -n r3.3 -c conda-forge mamba -y
#Installing Payload
# RUN mamba install -c conda-forge r-base=3.3.2

# RUN mamba install -c conda-forge r-plyr=1.8.3 r-foreach=1.4.3 r-doparallel=1.0.10 r-matrixstats=0.52.2 -y
# RUN mamba install -c bioconda bioconductor-deseq2=1.10.1 bioconductor-genomicranges=1.22.4 -y
# RUN R -e "tryCatch(install.packages('plyr', dependencies=TRUE, repos='https://cloud.r-project.org'),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(install.packages('doParallel', dependencies=TRUE, repos='https://cloud.r-project.org'),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(install.packages('foreach', dependencies=TRUE, repos='https://cloud.r-project.org'),  warning = function(w){ stop('install command gave a warning')})"

# # RUN R -e "install.packages('doParallel', version='1.0.10', dependencies=TRUE, repos='https://cloud.r-project.org')"  
# # RUN R -e "install.packages('foreach', version='1.4.3', dependencies=TRUE, repos='https://cloud.r-project.org')"
# RUN R -e "tryCatch(source('http://bioconductor.org/biocLite.R'),  warning = function(w){ stop('install command gave a warning')}); tryCatch(biocLite('GenomicRanges'),  warning = function(w){ stop('install command gave a warning')})" 
# RUN R -e "source('http://bioconductor.org/biocLite.R'); tryCatch(biocLite('DESeq2'),  warning = function(w){ stop('install command gave a warning')})" 



WORKDIR /MOUNT
#ENTRYPOINT SCRIPT
COPY --chown=nobody:nogroup ./ENTRYPOINT.sh /ENTRYPOINT.sh
COPY --chown=nobody:nogroup ./junctionseq.R /junctionseq.R
COPY --chown=nobody:nogroup ./dependencies.R /dependencies.R
RUN chmod 777 /ENTRYPOINT.sh
ENTRYPOINT ["/ENTRYPOINT.sh"]