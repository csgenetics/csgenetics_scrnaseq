################## BASE IMAGE #####################
FROM adoptopenjdk/openjdk11:ubuntu-jre-nightly 

################## METADATA #######################

LABEL base_image="adoptopenjdk/openjdk11:ubuntu-jre-nightly"
LABEL about.summary="Ubuntu image with Nextflow and Docker installed. Used to run the CS Genetics scRNA-Seq pipelin.e"

# Install required bash tools
RUN apt-get update && apt-get install -y apt-utils nano curl wget git awscli

# Install Nextflow
RUN cd /usr/local/bin && curl -s https://get.nextflow.io | bash

# Install Docker
RUN curl -fsSL https://get.docker.com -o get-docker.sh
RUN sh get-docker.sh

COPY entrypoint.sh /entrypoint.sh

ENTRYPOINT [ "/entrypoint.sh" ]


# Run nextflow 
# CMD ["/bin/bash"]
