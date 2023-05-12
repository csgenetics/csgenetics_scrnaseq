FROM python:3.9.7-slim

RUN apt-get update \
  && apt-get install -y libatlas-base-dev gfortran procps \
  && rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --upgrade pip

RUN pip install umi-tools-csgx==1.1.1
