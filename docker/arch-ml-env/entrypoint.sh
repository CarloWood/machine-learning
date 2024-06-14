#! /bin/bash

# Abort upon any error.
set -e

# Function to change ownership if the directory $1 is empty and not already owned by $2:$3.
# In that case it is initialize with the default files.
set_up_initial_machine_learning_dir() {
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
        # Also copy the initial env.source into its place.
        cp /usr/local/share/machine-learning/ml-env.source /home/archuser/machine-learning/env.source
        source /home/archuser/machine-learning/env.source
        mkdir "$GITACHE_ROOT"
      fi
    fi
  fi
}

# Fix ownership of /home/archuser/machine-learning and initialize it.
set_up_initial_machine_learning_dir /home/archuser/machine-learning $(id -u) $(id -g)

# Execute the original command.
exec "$@"
