# Build Docker image for summary_report module from this dockerfile 
FROM r-base

## install debian packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
libxml2-dev \
libcurl4-openssl-dev \
libssl-dev \
procps 

WORKDIR /script/

COPY install_packages.R /script/

RUN Rscript install_packages.R

CMD ["/bin/bash"]
