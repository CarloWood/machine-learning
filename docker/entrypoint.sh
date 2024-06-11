#! /bin/bash

# Abort upon any error.
set -e

# Update 
pacman --noconfirm -Syu
pacman-key --init

# Install sudo if it doesn't exist.
if ! command -v sudo &> /dev/null; then
  echo "Installing sudo..."
  pacman --noconfirm -Sy sudo
fi

# Create user archuser if it doesn't exist.
if ! id -u archuser > /dev/null 2>&1; then
  echo "Creating user archuser with uid $HOST_EUID and gid $HOST_EGID."
  groupadd --non-unique --gid $HOST_EGID archuser
  useradd --base-dir /home --no-create-home --gid $HOST_EGID --non-unique --uid $HOST_EUID -s /bin/bash archuser
fi

# Allow user archuser to use sudo without password.
if [ -d /etc/sudoers.d ]; then
  echo "archuser ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/archuser
else
  echo "archuser ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers
fi

# Make sure archuser's home directory is owned by archuser.
# This directory should exist as part of a docker volume; which is created by root.
chown archuser /home/archuser
chown archuser /home/archuser/data

# Switch to the archuser user and execute the original command.
exec su archuser -c "$@"
