################## BASE IMAGE #####################

FROM ubuntu:18.04

################## METADATA #######################

LABEL about.summary="Ubuntu image containing Python3 with HTML editing packages installed. Used to run the CS Genetics scRNA-Seq pipeline."

################## MAIN ###########################

RUN apt-get update && apt-get install -y apt-utils nano curl wget gosu 
RUN apt-get update && apt-get install -y software-properties-common gcc && \
    add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get update && apt-get install -y python3.6 python3-distutils python3-pip python3-apt


RUN pip3 install pandas \
jinja2 \
glob2 \
beautifulsoup4 

COPY entrypoint.sh /entrypoint.sh

# We run the entrypoint script to create a local user
# that has required docker permissions but a uid that
# matches that of the user who launched the docker container.
ENTRYPOINT [ "/entrypoint.sh" ]

# CMD will be replaced by the command line arguments.
CMD ["/bin/bash"]
