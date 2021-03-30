FROM bioconductor/bioconductor_docker:RELEASE_3_12

RUN R -e "devtools::install_github('JinmiaoChenLab/cytofkit', force=TRUE, depencencies=TRUE)"

ADD scripts/install_R_packages.R /tmp/
RUN Rscript /tmp/install_R_packages.R

ADD taproom/ /tmp/taproom/
RUN R -e "devtools::install('tmp/taproom')"
RUN R -e "install.packages('ggalluvial')"
#RUN R -e "devtools::install_github('camlab-bioml/taproom')"