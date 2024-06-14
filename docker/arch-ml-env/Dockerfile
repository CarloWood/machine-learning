FROM carlowood/arch-ml-bare:latest

RUN pacman --noconfirm -Sy which xterm

COPY skel /etc/skel
COPY cdeh /opt/cdeh

RUN useradd --create-home -s /bin/bash archuser

RUN echo "archuser ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers

RUN ln -s /opt/cdeh/addhistory /usr/local/bin/addhistory
RUN mkdir -p /usr/local/share/machine-learning
COPY ml-env.source /usr/local/share/machine-learning

# Set entrypoint script.
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Set the entrypoint.
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

# Change and set default user.
USER archuser

# Create /home/archuser/machine-learning as archuser.
RUN mkdir /home/archuser/machine-learning

# Make it the default working directory.
WORKDIR /home/archuser/machine-learning

# Declare volumes.
VOLUME ["/opt/cdeh/history", "/home/archuser/machine-learning", "/tmp/.X11-unix"]
