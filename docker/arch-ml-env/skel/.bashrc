# Don't put duplicate lines or lines starting with space in the history.
# See bash(1) for more options.
HISTCONTROL=ignoreboth

# Set a default environment.
export VERSION_CONTROL=numbered
export EDITOR=/usr/bin/nvim
export PATH="$HOME/bin:/usr/local/sbin:/usr/local/bin:/usr/bin:$HOME/.local/bin"
eval `dircolors --bourne-shell $HOME/.DIR_COLORS`

# Use 16 threads when using blas-openblas.
export OMP_NUM_THREADS=16

# Default top project for compiling libraries.
export TOPPROJECT=/home/archuser/machine-learning

# Aliases
alias vi="vim"
alias ls='ls --file-type --group-directories-first --color=auto'
# Print a list of env.source files.
alias pe='P=$(pwd); while : ; do E="$P/env.source"; if test -f "$E"; then echo $E; fi; [[ $P != "/" ]] || break; P=$(dirname "$P"); done'
# Print information on how to paste output to a paste site from the command line.
alias ix='echo "<command to print output> |& curl -F '"'"'f:1=<-'"'"' ix.io"'

# Set a fancy color prompt
PS1PREFIX="docker:"
PS1=$PS1PREFIX'\[\e[35m\]\w\[\e[32`if [ $EUID = 0 ]; then echo ";7"; fi`m\]>\[\e[30;47;0m\]'

# Check the window size after each command and, if necessary,
# update the values of LINES and COLUMNS.
shopt -s checkwinsize

# Show aliases and bash functions when invoking `which`.
function which()
{
  (alias; declare -f) | /usr/bin/which --tty-only --read-alias --read-functions --show-tilde --show-dot "$@"
}

# This must be done as the VERY last thing in .bashrc!
#export CDEH_VERBOSE=1
CDEH_ROOT=/opt/cdeh
. $CDEH_ROOT/env.bashrc
