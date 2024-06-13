# This DockerFile prepares an archlinux container with only all desired packages already installed.
# Lets call it arch-ml-bare.
#
FROM amd64/archlinux:latest

RUN pacman --noconfirm -Syu \
    && pacman-key --init \
    && pacman --noconfirm -Sy sudo git clang make pkgconfig cmake boost pango cairo eigen tensorflow-opt ttf-dejavu neovim ccache
