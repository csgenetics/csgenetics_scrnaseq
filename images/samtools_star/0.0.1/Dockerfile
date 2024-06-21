FROM continuumio/miniconda3
# https://medium.com/@chadlagore/conda-environments-with-docker-82cdc9d25754

RUN apt-get update && apt-get install unzip build-essential libz-dev libncurses5-dev libncursesw5-dev liblzma-dev libbz2-dev libcurl4-openssl-dev -y
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.18/samtools-1.18.tar.bz2/download -O samtools-1.18.tar.bz2
RUN tar xvjf samtools-1.18.tar.bz2
RUN cd samtools-1.18 && make
ENV PATH $PATH:/samtools-1.18

ADD star_env.yml /tmp/star_env.yml
RUN conda env create -f /tmp/star_env.yml
# This will not have any effect for Singularity containers
RUN echo "source activate star_env" > ~/.bashrc
# Rather, manually append the env to the path so it is the default python run
ENV PATH /opt/conda/envs/star_env/bin:$PATH