#!/bin/bash

set -e

# Use specified LOCAL_USER_ID to create a new user
# that matches the USER ID of the user launching the container
# See https://denibertovic.com/posts/handling-permissions-with-docker-volumes/
# for further details on motives to work with the local user's UID.
USER_ID=${LOCAL_USER_ID:-9001}

# We also need the gid of the docker group so that we can
# add the user to that grop so that they still have access to
# the docker.sock file in
DOCKER_GID=${DOCKER_GID}

echo "Starting with UID : $USER_ID"

# Create the user's group
addgroup --gid $USER_ID user

# Create the user and associate to the group
useradd --shell /bin/bash -u $USER_ID -o -c "" -g $USER_ID -m user

# get the docker group id that we need to add the user to
DOCKER_GID=$(ls -l /var/run/docker.sock | cut -d" " -f4)

# Try to add the user to that group
# If group doesn't exist, create the group and call it docker_sock
# Then add the user to that group
usermod -aG $DOCKER_GID user || (addgroup --gid $DOCKER_GID docker_sock; echo "Creating group docker_sock with gid ${DOCKER_GID}"; usermod -aG $DOCKER_GID user; echo "Added user to group docker_sock")

# usermod -aG docker user
export HOME=/home/user

# exec and run the actual process specified in the
# CMD of the Dockerfile (which gets passed as ${@})
# using gosu so that it is run as the user.
# The user should have the required groups to run the /var/run/docker.sock
# so that subsequent docker containers can be created
exec gosu user "$@"
