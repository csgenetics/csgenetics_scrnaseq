## Build a image from this dockerfile for having nextflow and docker already installed in container used it for dockerized running version
FROM adoptopenjdk/openjdk11:ubuntu-jre-nightly
# Docker Image having Java tool as base image along with nextflow and docker installed inside the container

# Install required bash tools
RUN apt-get update && apt-get install -y apt-utils nano curl wget git awscli

# Install Nextflow
RUN cd /usr/local/bin && curl -s https://get.nextflow.io | bash

# Install Docker
RUN curl -fsSL https://get.docker.com -o get-docker.sh
RUN sh get-docker.sh
CMD ["/bin/bash", "nextflow"]
