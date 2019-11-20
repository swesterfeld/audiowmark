#!/usr/bin/env python3

# test how long the watermarker takes until the first audio sample is available

import subprocess
import time

for i in range (10):
    start_time = time.time()
    proc = subprocess.Popen (["audiowmark", "add", "test/T01__09_sing_sang_sung.flac.wav", "-", "f0"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # we wait for actual audio data, so we read somewhat larger amount of data that the wave header
    x = proc.stdout.read (1000)
    end_time = time.time()
    print (end_time - start_time, x[0:4], len (x))
