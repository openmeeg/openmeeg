#
# find openmeeg lib
#
import sys
import os

# import locally built openmeeg
topdir = os.getcwd()
omdir = ""

# heuristic
found = False
l = ""
for i in range(1, 10):
    d = os.path.join( topdir, l , "build" , "wrapping", "src")
    if os.path.exists( os.path.join( d , "_openmeeg.so")):
        omdir = os.path.join( topdir, l)
        sys.path.append(d)
        found = True
        break
    l = os.path.join(l , "..")

if not found:
    print("Cannot find locally compiled _openmeeg.so file")
    print("Please goto openmeeg repository or a subdir of it")
    sys.exit(-1)
