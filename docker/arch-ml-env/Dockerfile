FROM carlowood/arch-ml-bare:latest

RUN pacman --noconfirm -Sy which

COPY skel /etc/skel
COPY cdeh /opt/cdeh

RUN useradd --create-home -s /bin/bash archuser

RUN echo "archuser ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers

RUN ln -s /opt/cdeh/addhistory /usr/local/bin/addhistory

# Change and set default user.
USER archuser

# Create /home/archuser/machine-learning as archuser.
RUN mkdir /home/archuser/machine-learning
# Put the environment of the project in there.
COPY ml-env.source /home/archuser/machine-learning

# Make it the default working directory.
WORKDIR /home/archuser/machine-learning
