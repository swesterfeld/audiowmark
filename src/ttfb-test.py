#!/usr/bin/env python3

# test how long the watermarker takes until the first audio sample is available

import subprocess
import shlex
import time
import sys

seconds = 0

for i in range (10):
    start_time = time.time() * 1000
    proc = subprocess.Popen (shlex.split (sys.argv[1]), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # we wait for actual audio data, so we read somewhat larger amount of data that the wave header
    x = proc.stdout.read (1000)
    end_time = time.time() * 1000

    seconds += end_time - start_time
    print ("%.2f" % (end_time - start_time), x[0:4], len (x))

print ("%.2f" % (seconds / 10), "avg")
