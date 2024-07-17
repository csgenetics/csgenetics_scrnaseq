#!/bin/bash

## Define GCP parameters
GCP_PROJECT="eric-sandbox-421120"
GCP_REGION="us-central1"
GCP_REGISTRY="us-central1-docker.pkg.dev/eric-sandbox-421120/csgenetics"

## Define image names
RUN_NAME="nextflow-run-test" # Must be lowercase, no underscores
IMAGE_NAME="csgenetics-nextflow-run:${RUN_NAME}"
GCP_IMAGE="${GCP_REGISTRY}/${IMAGE_NAME}"

## Build image locally
docker build --no-cache -t ${IMAGE_NAME} -f images/gcp_launch/Dockerfile .

## Tag and push to GCP
docker tag ${IMAGE_NAME} ${GCP_IMAGE}
docker push ${GCP_IMAGE}

## Run batch job
echo "Running batch job..."
gcloud beta batch jobs submit ${RUN_NAME} --project ${GCP_PROJECT} --location ${GCP_REGION} --config - <<EOD
{
  "name": "projects/${GCP_PROJECT}/locations/${GCP_REGION}/jobs/${RUN_NAME}",
  "taskGroups": [
    {
      "taskCount": "1",
      "parallelism": "1",
      "taskSpec": {
        "computeResource": {
          "cpuMilli": "1000",
          "memoryMib": "512"
        },
        "runnables": [
          {
            "container": {
              "imageUri": "${GCP_IMAGE}",
              "entrypoint": "",
              "volumes": []
            }
          }
        ],
        "volumes": []
      }
    }
  ],
  "allocationPolicy": {
    "instances": [
      {
        "policy": {
          "provisioningModel": "STANDARD",
          "machineType": "e2-micro"
        }
      }
    ]
  },
  "logsPolicy": {
    "destination": "CLOUD_LOGGING"
  }
}
EOD