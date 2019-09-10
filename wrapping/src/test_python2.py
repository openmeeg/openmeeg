#!/usr/bin/env python

#
import sys
sys.path.append("/Users/jls/Development/athena/openmeeg/build/wrapping/src")


import openmeeg as om


dipoles = om.Matrix()
dipoles.load("/Users/jls/Development/athena/openmeeg/data/Head1/Head1.dip")
D = om.asarray(dipoles)
