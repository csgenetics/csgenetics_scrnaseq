################## BASE IMAGE #####################

FROM ubuntu:22.04

################## METADATA #######################

LABEL about.summary="Ubuntu image containing the AWS CLI. Used for downloading files from AWS S3 buckets."

################## MAIN ###########################
RUN apt-get update && apt-get install -y curl unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

# CMD will be replaced by the command line arguments.
CMD ["/bin/bash"]