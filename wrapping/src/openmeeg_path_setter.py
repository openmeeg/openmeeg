#
# find openmeeg lib
#
import sys
import os

# import locally built openmeeg
topdir = os.getcwd()
omdir  = os.path.join(topdir, "build/wrapping/src/")
if not os.path.exists( os.path.join(omdir , "_openmeeg.so")):
    print("Go to openmeeg topdir before launching this script")
    exit(-1)

sys.path.append(omdir)
