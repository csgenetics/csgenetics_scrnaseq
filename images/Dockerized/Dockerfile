################## BASE IMAGE #####################
FROM adoptopenjdk/openjdk11:ubuntu-jre-nightly 

################## METADATA #######################

LABEL base_image="adoptopenjdk/openjdk11:ubuntu-jre-nightly"
LABEL about.summary="Ubuntu image with Nextflow and Docker installed. Used to run the CS Genetics scRNA-Seq pipeline."

# Install required bash tools
RUN apt-get update && apt-get install -y apt-utils nano curl wget git awscli gosu

# Install Nextflow
RUN cd /usr/local/bin && curl -s https://get.nextflow.io | bash
RUN chmod o+rx /usr/local/bin/nextflow

# Install Docker
RUN curl -fsSL https://get.docker.com -o get-docker.sh
RUN sh get-docker.sh

COPY entrypoint.sh /entrypoint.sh

# We run the entrypoint script to create a local user
# that has required docker permissions but a uid that
# matches that of the user who launched the docker container.
ENTRYPOINT [ "/entrypoint.sh" ]

# CMD will be replaced by the command line arguments.
CMD ["/bin/bash"]
