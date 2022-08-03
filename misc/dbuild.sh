#!/bin/bash
set -Eeuo pipefail

docker build -f "misc/Dockerfile" -t audiowmark-dbuild .
docker build -f "misc/Dockerfile-arch" -t audiowmark-dbuild-arch .
