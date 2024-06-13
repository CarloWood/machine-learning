#! /bin/bash

# Abort upon any error.
set -e

# Make sure archuser's history directory is owned by archuser.
# This directory should exist as part of a docker volume; which is created by root.
sudo chown -R archuser:archuser /opt/cdeh/history/archuser

# Fix ownership of home directory.
sudo chown -R archuser:archuser /home/archuser

# Execute command.
tail -f /dev/null
