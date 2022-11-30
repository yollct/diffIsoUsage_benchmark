FROM rocker/r-ver:3.3


RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools unar libbz2-dev liblzma-dev libglpk-dev git 
RUN apt update && apt install parallel -y  --force-yes
# Installing Mamba
RUN mkdir /Rlib
RUN chmod 777 /Rlib

RUN echo "R_LIBS_USER='/Rlib'" >> /usr/local/lib/R/etc/Renviron
RUN echo "R_LIBS='/Rlib'" >> /usr/local/lib/R/etc/Renviron
# RUN R -e "tryCatch(install.packages(\"tidyverse\", lib=\"/Rlib"),  warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(source(\"http://bioconductor.org/biocLite.R\"),  warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(BiocManager::install(\"DESeq2\",ask=FALSE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(biocLite(\"JunctionSeq\"), warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(biocLite(\"DESeq2\"), warning = function(w){ stop(\"install command gave a warning\")})"

WORKDIR /MOUNT
#ENTRYPOINT SCRIPT
COPY --chown=nobody:nogroup ./ENTRYPOINT.sh /ENTRYPOINT.sh
COPY --chown=nobody:nogroup ./junctionseq.R /junctionseq.R
COPY --chown=nobody:nogroup ./dependencies.R /dependencies.R
RUN chmod 777 /ENTRYPOINT.sh
ENTRYPOINT ["/ENTRYPOINT.sh"]
