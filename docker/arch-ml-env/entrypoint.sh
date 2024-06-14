#! /bin/bash

# Abort upon any error.
set -e

# Function to change ownership if the directory $1 is empty and not already owned by $2:$3.
fix_ownership_if_empty_and_not_owned() {
  local dir=$1
  local uid=$2
  local gid=$3

  # Check if the directory exists.
  if [ -d "$dir" ]; then
    # Get current owner and group of the directory.
    current_uid=$(stat -c "%u" "$dir")
    current_gid=$(stat -c "%g" "$dir")

    # Check if the directory is empty and not owned by the specified uid and gid.
    if [ "$current_uid" != "$uid" ] || [ "$current_gid" != "$gid" ]; then
      if [ ! "$(ls -A $dir)" ]; then
        echo "Directory $dir is empty and not owned, changing ownership to $uid:$gid."
        sudo chown $uid:$gid $dir
      fi
    fi
  fi
}

# Fix ownership of /home/archuser/machine-learning.
fix_ownership_if_empty_and_not_owned /home/archuser/machine-learning $(id -u) $(id -g)

# Execute the original command.
exec "$@"
