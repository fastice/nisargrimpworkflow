#!/bin/bash
#
# Cron entry point for autoupdateNISAR.
#
# cron runs with an almost-empty environment, so this wrapper sets up the two
# things the daily update needs on PATH -- the conda env (autoupdateNISAR,
# searchASF, ariaDownload, SetupNISAR, ... plus GDAL libs) and the GrIMP C
# binaries (simoffsets, cullst, siminsar, mosaic3d) -- then runs the update for
# one project directory. It lives in the package so it is versioned with the code;
# reference it from crontab by its absolute path.
#
# Usage (from cron, by absolute path):
#   runAutoupdateNISAR.sh <projectDir> [extra autoupdateNISAR args]
# e.g.
#   runAutoupdateNISAR.sh /Volumes/insar4/ian/Antarctica2026/Antarctica40
#   runAutoupdateNISAR.sh /Volumes/.../Antarctica40 --noDownload
#
# Environment overrides (defaults suit the GrIMP workstation):
#   CONDA_BASE  conda/miniforge install root    [/home/ian/miniforge3]
#   BINDIR      dir holding the GrIMP C binaries [/home/ian/bin/x86_64]
#
# Note: no `set -u` here -- conda's own activate/deactivate hooks reference
# unbound variables, so `set -u` would abort activation. The explicit guards
# below cover the missing-argument and unmounted-tree cases instead.

CONDA_BASE="${CONDA_BASE:-/home/ian/miniforge3}"
BINDIR="${BINDIR:-/home/ian/bin/x86_64}"

if [ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    echo "$(date): cannot find conda at $CONDA_BASE (set CONDA_BASE)" >&2
    exit 1
fi
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate base
export PATH="$BINDIR:$PATH"

PROJ="${1:-}"
if [ -z "$PROJ" ]; then
    echo "usage: $0 <projectDir> [extra autoupdateNISAR args]" >&2
    exit 2
fi
shift

# Bail cleanly if the (NFS) project tree is not mounted / ready, so a failed
# mount does not turn into a spurious "no data" run.
if [ ! -f "$PROJ/autoupdate.yaml" ]; then
    echo "$(date): no autoupdate.yaml in $PROJ (mount not ready?)" >&2
    exit 1
fi

cd "$PROJ" || exit 1
# The crontab line redirects into dailyDownloadLogs/cron.log; the shell opens
# that redirect BEFORE anything runs, so the directory must exist up front on
# a fresh project. Harmless if already present.
mkdir -p dailyDownloadLogs
exec autoupdateNISAR "$@" autoupdate.yaml
