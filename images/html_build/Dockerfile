FROM continuumio/miniconda3

ADD html_build.yml /tmp/html_build.yml
RUN conda env create -f /tmp/html_build.yml

# This will not have any effect for Singularity containers
RUN echo "source activate html_build" > ~/.bashrc
# Rather, manually append the env to the path so it is the default python run
ENV PATH /opt/conda/envs/html_build/bin:$PATH