#!/bin/bash

set -e

# use specified user name or use `default` if not specified
MY_USERNAME="${MY_USERNAME:-default}"
echo "MY_USERNAME is set to ${MY_USERNAME}"

# use specified group name or use the same user name also as the group name
MY_GROUP="${MY_GROUP:-${MY_USERNAME}}"
echo "MY_GROUP is set to ${MY_GROUP}"

# use the specified UID for the user
MY_UID="${MY_UID:-1000}"
echo "MY_UID is set to ${MY_UID}"

# use the specified GID for the user
MY_GID="${MY_GID:-${MY_UID}}"
echo "MY_GID is set to ${MY_GID}"

# check to see if group exists; if not, create it
if grep -q -E "^${MY_GROUP}:" /etc/group > /dev/null 2>&1
then
  echo "INFO: Group exists; skipping creation"
else
  echo "INFO: Group doesn't exist; creating..."
  # create the group
  addgroup --gid $MY_GID $MY_GROUP || (echo "INFO: Group exists but with a different name; renaming..."; groupmod -g "${MY_GID}" -n "${MY_GROUP}" "$(awk -F ':' '{print $1":"$3}' < /etc/group | grep ":${MY_GID}$" | awk -F ":" '{print $1}')")
fi


# check to see if user exists; if not, create it
if id -u "${MY_USERNAME}" > /dev/null 2>&1
then
  echo "INFO: User exists; skipping creation"
else
  echo "INFO: User doesn't exist; creating..."
  # create the user
  adduser --uid $MY_UID --gid $MY_GID --home "/home/${MY_USERNAME}" --shell /bin/bash --disabled-password --gecos "" $MY_USERNAME
fi

# exec and run the actual process specified in the CMD of the Dockerfile (which gets passed as ${@})
exec "$@"

