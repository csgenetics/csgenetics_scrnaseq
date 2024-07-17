# Build the image for running nextflow
FROM nextflow/nextflow:24.04.3

## Use root & switch context to /app directory
USER root
WORKDIR /app

## Copy the modified repo from the current working directory
COPY ./ /app/csgenetics_scrnaseq

## Run the nextflow pipeline
CMD ["nextflow", "run", "/app/csgenetics_scrnaseq/main.nf", "-profile", "gcp_test", "-resume"]