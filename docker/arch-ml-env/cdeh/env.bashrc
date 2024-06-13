# cdeh v2.0 -- an automatic history and environment switch system.
#
# Copyright (C) 2007 (v1.0), 2018 (v2.0)  Carlo Wood <carlo@alinoe.com>

# Write history on exit.
function exit ()
{
  history -a $HISTFILE
  unset -f exit
  rm -rf "/tmp/cdeh.$CDEH_USER/$$"
  exit $1
}

function __cdeh_sanity ()
{
  # If CDEH_ROOT is not set, print a warning and do nothing.
  if test -z "$CDEH_ROOT"; then
    echo "cdeh: CDEH_ROOT is not set prior to sourcing env.bashrc. cdeh disabled." >&2
    return 0
  fi
  if test ! -d "$CDEH_ROOT"; then
    echo "cdeh: CDEH_ROOT = \"$CDEH_ROOT\"" 2>&1
    echo "cdeh: $CDEH_ROOT: No such directory. cdeh disabled." 2>&1
    return 0
  fi
  if test ! -d "$CDEH_ROOT/history"; then
    echo "cdeh: CDEH_ROOT = \"$CDEH_ROOT\"" 2>&1
    echo "cdeh: $CDEH_ROOT/history: No such directory. cdeh disabled." 2>&1
    return 0
  fi
  if test ! -d "$CDEH_ROOT/history/$CDEH_USER"; then
    if ! mkdir "$CDEH_ROOT/history/$CDEH_USER"; then
      echo "cdeh: Disabled. Please run as root: chmod 1777 $CDEH_ROOT/history" 2>&1
      return 0
    fi
  fi
  return 1
}

export CDEH_USER="`/usr/bin/whoami`"
# Do some sanity checks.
if __cdeh_sanity; then
  unset -f exit
  # This is just to clear PROMPT_COMMAND in case an error was
  # introduced after first having succesfully set PROMPT_COMMAND.
  CDEH_TEST=`echo $PROMPT_COMMAND | grep do_prompt`;
  if test -n "$CDEH_TEST"; then
    unset PROMPT_COMMAND
  fi
  unset CDEH_TEST
  unset CDEH_USER
else

# Keep backup of the value of HISTSIZE if any.
# Default value is 1000.
export CDEH_HISTSIZE=${HISTSIZE-1000}
# Erase history
HISTFILE=
HISTSIZE=0
# Set sane value for HISTFILESIZE
test ${HISTFILESIZE-0} -lt $CDEH_HISTSIZE && HISTFILESIZE=$CDEH_HISTSIZE
export HISTFILESIZE
export CDEH_TMP="/tmp/cdeh.$CDEH_USER/$$"
export CDEH_HISTROOT="$CDEH_ROOT/history/$CDEH_USER"
mkdir -p "$CDEH_TMP"
echo -n > "$CDEH_TMP/prevwd"
echo -n > "$CDEH_TMP/prevhd"
echo -n > "$CDEH_TMP/preved"

function __cdeh_store_environment ()
{
  declare -p | awk '
      BEGIN { show=0 }
      /^declare -/ {
        name=$3;
        gsub("=.*","",name);
        gsub("-+", "-g", $2);
        show=!index($2,"r") && !match(name, "^((__cdeh_|CDEH_|BASH|PULSE_PROP_OVERRIDE)|(name|HISTFILE|HISTFILESIZE|HISTSIZE|PIPESTATUS|PWD|_|SSH_AGENT_PID|SSH_AUTH_SOCK)$)")
      }
      {
        if (show)
          print
      }' > "$CDEH_TMP/env.base"
  declare -pf | awk '
      BEGIN { show=0 }
      /^[[:alnum:]_-]+ \(\)/ {
        name=$1;
        show=!match(name, "^__cdeh_")
      }
      /^declare -/ {
        name=$3;
        gsub("-+", "-g", $2);
        show=!match(name, "^__cdeh_") }
      {
        if (show)
          print
      }' >> "$CDEH_TMP/env.base"
  alias -p >> "$CDEH_TMP/env.base"
}

function __cdeh_clear_environment ()
{
  test -z "$CDEH_VERBOSE" || echo "$$ CLEARING ENVIRONMENT"
  for name in $(declare -p | \
                /bin/grep '^declare -[AFafgilntux-]* ' | \
                /bin/sed -e 's/=.*//;s/.* //' | \
                /bin/grep -E -v '^((__cdeh_|CDEH_|BASH_|PULSE_PROP_OVERRIDE)|(BASH|BASHPID|BASHOPTS|COMP_WORDBREAKS|DIRSTACK|FUNCNAME|GROUPS|LINENO|RANDOM|SECONDS|HOME|PATH|PS1|PROMPT_COMMAND|PWD|name|HISTFILE|HISTFILESIZE|HISTSIZE|HISTCMD|HISTCONTROL|PIPESTATUS|HOSTFILE|MAILCHECK|_|SSH_AGENT_PID|SSH_AUTH_SOCK)$)'); do
    unset $name
  done
  for name in $(declare -F | \
                /bin/grep '^declare -[[:alnum:]]*f' | \
                /bin/sed -e 's/.* //' | \
                /bin/grep -E -v '^__cdeh_'); do
    unset -f $name
  done
  for name in $(alias -p | \
                /bin/grep '^alias [[:alnum:]_-]*=' | \
                /bin/sed -e 's/=.*//;s/^alias //'); do
    unalias $name
  done
  unset name
}

function __cdeh_resource ()
{
  CDEH_SAVE_IFS="$IFS"
  IFS=$'\n'
  for __cdeh_env_file in $CDEH_ENVFILES; do
    if test "$__cdeh_env_file" = "__cdeh_reset_environment"; then
      __cdeh_clear_environment
      test -z "$CDEH_VERBOSE" || echo "$$ RELOADING SAVED ENVIRONMENT from $CDEH_TMP/env.base"
      source "$CDEH_TMP/env.base"
    else
      test -z "$CDEH_VERBOSE" || echo "$$ Sourcing \"$__cdeh_env_file\""
      source "$__cdeh_env_file"
    fi
  done
  IFS="$CDEH_SAVE_IFS"
  unset CDEH_SAVE_IFS
}

function __cdeh_rebuild_source_dirs ()
{
  CDEH_CURRENT_DIR="/"
  cdeh_i=0
  unset CDEH_SOURCE_DIRS
  CDEH_DESCENT_DIR="$PWD"
  until [ "$CDEH_DESCENT_DIR" = "" ]; do
    if [ -f "$CDEH_DESCENT_DIR/env.source" ]; then
      CDEH_SOURCE_DIRS[cdeh_i]="$CDEH_DESCENT_DIR"
      cdeh_i=$((cdeh_i + 1))
    fi
    CDEH_DESCENT_DIR="${CDEH_DESCENT_DIR%/*}"
  done
  if test $cdeh_i -gt 0; then
    CDEH_CURRENT_DIR="${CDEH_SOURCE_DIRS[0]}"
    CDEH_NEW_DIR="$CDEH_CURRENT_DIR"
  else
    CDEH_NEW_DIR="/"
  fi
  if [ -f "/env.source" ]; then
    CDEH_SOURCE_DIRS[cdeh_i]=""
    cdeh_i=$((cdeh_i + 1))
    if test $cdeh_i -eq 1; then
      CDEH_NEW_DIR="/"
    fi
  fi
}

function __cdeh_rebuild_envfiles ()
{
  echo "__cdeh_reset_environment" >> "$CDEH_TMP/envfiles"
  cdeh_j=$cdeh_i
  while test $cdeh_j -gt 0; do
   cdeh_j=$((cdeh_j - 1))
   echo ${CDEH_SOURCE_DIRS[$cdeh_j]}/env.source >> "$CDEH_TMP/envfiles"
  done
}

declare -xf __cdeh_store_environment __cdeh_resource __cdeh_rebuild_source_dirs __cdeh_rebuild_envfiles

function resource ()
{
  __cdeh_rebuild_source_dirs
  rm -f "$CDEH_TMP/envfiles"
  __cdeh_rebuild_envfiles
  export CDEH_ENVFILES=$(cat "$CDEH_TMP/envfiles")
  __cdeh_resource
}

function no_path()
{
  eval "case :\$${2-PATH}: in *:$1:*) return 1;; *) return 0;; esac"
}

function add_path ()
{
  [ -d ${1:-.} ] && no_path $* && eval ${2:-PATH}="\$${2:-PATH}:$1"
}

function pre_path ()
{
  [ -d ${1:-.} ] && no_path $* && eval ${2:-PATH}="$1:\$${2:-PATH}"
}

function del_path ()
{
  no_path $* || eval ${2:-PATH}=`eval echo :'$'${2:-PATH}: | sed -e "s;:$1:;:;g" -e "s;^:;;" -e "s;:\$;;"`
}

export TOPPROJECT=${TOPPROJECT=TOPPROJECT_IS_NOT_SET}

if [ "$TOPPROJECT" = "TOPPROJECT_IS_NOT_SET" -a ! -z "$DEFAULT_TOPPROJECT" ]; then
  echo "Setting TOPPROJECT to \"$DEFAULT_TOPPROJECT\"!"
  export TOPPROJECT=$DEFAULT_TOPPROJECT
fi

export PROMPT_COMMAND='
    if "$CDEH_ROOT"/do_prompt $$ "$HISTFILE"; then
      if test -f "$CDEH_TMP/histfile"; then
        HISTSIZE=0;
        HISTSIZE=$CDEH_HISTSIZE;
        export HISTFILE=$(/bin/cat "$CDEH_TMP/histfile");
        history -r;
      fi;
      if test -f "$CDEH_TMP/envfiles"; then
        export CDEH_ENVFILES=$(/bin/cat "$CDEH_TMP/envfiles");
        __cdeh_resource;
      fi;
    fi'

# Store the environment as it is at the end of ~/.bashrc.
export CDEH_STORE_ENVIRONMENT="yes"

fi # __cdeh_sanity succeeded.
