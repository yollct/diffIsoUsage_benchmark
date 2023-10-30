FROM rocker/r-ver:4.2

RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools unar libbz2-dev liblzma-dev libglpk-dev libnetcdf-c++4-1 git libfontconfig1-dev cmake python3
RUN apt update && apt install parallel -y  --force-yes
# Installing Mamba
RUN mkdir /Rlib
RUN chmod 777 /Rlib
COPY ./rSeqDiff.beta.0.1.tar.gz /rSeqDiff.beta.0.1.tar.gz
COPY ./rseq-0.2.2-src /rseq-0.2.2-src
COPY ./code.R /code.R

RUN echo "R_LIBS_USER='/Rlib'" >> /usr/local/lib/R/etc/Renviron
RUN echo "R_LIBS='/Rlib'" >> /usr/local/lib/R/etc/Renviron
RUN R -e "tryCatch(install.packages(\"BiocManager\", repo=\"https://cloud.r-project.org\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"dplyr\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"tidyr\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"ggplot2\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"Rcpp\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(install.packages(\"/ncdf4_1.21.tar.gz\", repos=NULL, type=\"source\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"

#RUN R -e "tryCatch(install.packages(\"car\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(BiocManager::install(\"Rhdf5lib\",ask=FALSE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(BiocManager::install(\"tximport\",ask=FALSE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(BiocManager::install(\"BiocParallel\",ask=FALSE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(BiocManager::install(\"mzR\",ask=FALSE, dependencies=TRUE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(BiocManager::install(\"BiocGenerics\",ask=FALSE, dependencies=TRUE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(BiocManager::install(\"ProtGenerics\",ask=FALSE, dependencies=TRUE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
# RUN R -e "tryCatch(BiocManager::install(\"Biobase\",ask=FALSE, dependencies=TRUE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"

#RUN R -e "tryCatch(install.packages(\"/mzR_2.34.1.tar.gz\", repos=NULL, type=\"source\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"/rSeqDiff.beta.0.1.tar.gz\", repos=NULL, type=\"source\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(BiocManager::install(\"Rsamtools\", version=\"devel\", ask=FALSE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"

#RUN R -e "tryCatch(install.packages(\"/NBSplice_1.15.0.tar.gz\", repos=NULL, type=\"source\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(biocLite(\"JunctionSeq\"), warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(biocLite(\"DESeq2\"), warning = function(w){ stop(\"install command gave a warning\")})"

WORKDIR /MOUNT
#ENTRYPOINT SCRIPT
COPY --chown=nobody:nogroup ./ENTRYPOINT.sh /ENTRYPOINT.sh
COPY --chown=nobody:nogroup ./rseqdiff.R /rseqdiff.R

RUN chmod 777 /ENTRYPOINT.sh
ENTRYPOINT ["/ENTRYPOINT.sh"]
