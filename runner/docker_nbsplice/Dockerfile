FROM rocker/r-ver:4.2


RUN apt-get update && apt-get --no-install-recommends --fix-broken install -y curl libcurl4-openssl-dev libssl-dev libxml2-dev vim samtools unar libbz2-dev liblzma-dev libglpk-dev git libfontconfig1-dev cmake
RUN apt update && apt install parallel -y  --force-yes
# Installing Mamba
RUN mkdir /Rlib
RUN chmod 777 /Rlib
COPY ./NBSplice_1.15.0.tar.gz /NBSplice_1.15.0.tar.gz
COPY ./mppa_1.0.tar.gz /mppa_1.0.tar.gz

RUN echo "R_LIBS_USER='/Rlib'" >> /usr/local/lib/R/etc/Renviron
RUN echo "R_LIBS='/Rlib'" >> /usr/local/lib/R/etc/Renviron
RUN R -e y"trCatch(install.packages(\"BiocManager\", repo=\"https://cloud.r-project.org\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"dplyr\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"tidyr\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"ggplot2\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"reshape2\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"car\", repo=\"https://cloud.r-project.org\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"tximport\",ask=FALSE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"BiocParallel\",ask=FALSE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(BiocManager::install(\"edgeR\",ask=FALSE, dependencies=TRUE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"/mppa_1.0.tar.gz\", repos=NULL, type=\"source\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(BiocManager::install(\"NBSplice\", version=\"devel\", ask=FALSE, force=TRUE, lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
RUN R -e "tryCatch(install.packages(\"/NBSplice_1.15.0.tar.gz\", repos=NULL, type=\"source\", lib=\"/Rlib\"),  warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(biocLite(\"JunctionSeq\"), warning = function(w){ stop(\"install command gave a warning\")})"
#RUN R -e "tryCatch(biocLite(\"DESeq2\"), warning = function(w){ stop(\"install command gave a warning\")})"

WORKDIR /MOUNT
#ENTRYPOINT SCRIPT
COPY --chown=nobody:nogroup ./ENTRYPOINT.sh /ENTRYPOINT.sh
COPY --chown=nobody:nogroup ./NBsplice.R /NBsplice.R
COPY --chown=nobody:nogroup ./dependencies.R /dependencies.R
RUN chmod 777 /ENTRYPOINT.sh
ENTRYPOINT ["/ENTRYPOINT.sh"]
