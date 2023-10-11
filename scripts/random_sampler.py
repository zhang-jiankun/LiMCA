#!/bin/bash

import sys
import random

def random_sampler(fn, k):
    """
        See http://metadatascience.com/2014/02/27/random-sampling-from-very-large-files/
    """
    sample = []
    k = int(k)
    random.seed(100)
    with open(fn, 'rb') as f:
		f.seek(0, 2)
		filesize = f.tell()
		random_set = sorted(random.sample(xrange(filesize), k))
		for i in xrange(k):
			f.seek(random_set[i])
			# Skip current line (because we might be in the middle of a line) 
			f.readline()
			# Append the next line to the sample set 
			sample.append(f.readline().rstrip())
    return sample

sample = random_sampler(*sys.argv[1:])
for line in sample:
    sys.stdout.write(line + '\n')
