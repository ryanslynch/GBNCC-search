#!/bin/env sh

export OUTDIR=`awk '/^topdir/{print $3}' pipeline/lib/python/config.py | tr -d \"`

download.py > $OUTDIR/download.out &
submit.py > $OUTDIR/submit.out &
trackjobs.py > $OUTDIR/trackjobs.out &
