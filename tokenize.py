#! /usr/bin/env python

import sys

n = len(sys.argv)
print("token count: %d" % n)
for i in range(n):
    print("token[%d]: %s" %(i, sys.argv[i]))
