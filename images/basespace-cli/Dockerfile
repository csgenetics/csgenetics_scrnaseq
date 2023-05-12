FROM ubuntu:bionic

RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

RUN wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /bin/bs

RUN chmod u+x /bin/bs
