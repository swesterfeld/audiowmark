#!/bin/bash

# Example script how to replace the audio stream in a MOV
# container with watermarked audio
#
# This file is in the public domain
#
# Author: Nils Goroll <nils.goroll@uplex.de>

function usage {
    cat >&2 <<EOF
Usage: $0 <watermark> <inputfile> <outputfile>
EOF
}

typeset -A cmd
cmd[audiowmark]=$(which audiowmark)
cmd[ffmpeg]=$(which ffmpeg)

if [[ "${#@}" -ne "3" ]] ; then
    usage
    exit 1
fi

if ! [[ -f "${2}" ]] ; then
    echo >&2 Input "${2}" is not a file
    exit 1
fi

for c in "${!cmd[@]}" ; do
    if [[ -z "${cmd[$c]}" ]] ; then
	echo >&2 Required tool "${c}" not found
	exit 1
    fi
done

set -eu

atmp=$(mktemp).mov
mtmp=$(mktemp).wav
vtmp=$(mktemp).mov

# split audio/video
${cmd[ffmpeg]} -i "${2}" \
	-map 0:0 -vcodec copy ${vtmp} \
	-map 0:1 -acodec copy ${atmp}

# watermark audio it into ${mtmp}
${cmd[ffmpeg]} -i ${atmp} -ac 2 -f wav pipe:1 |
    ${cmd[audiowmark]} add - ${mtmp} "${1}" || true

# join into mov container
${cmd[ffmpeg]} -i ${vtmp} -i ${mtmp} -c:v copy -c:a aac "${3}"
rm -f ${atmp} ${mtmp} ${vtmp}
