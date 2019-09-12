#
# find openmeeg lib
#
import sys
import os

# import locally built openmeeg
topdir = os.getcwd()

# heuristic
found = 0;
l = ""
for i in range(1, 10):
    d = os.path.join( topdir, l , "build" , "wrapping", "src")
    if os.path.exists( os.path.join( d , "_openmeeg.so")):
        print(d)
        sys.path.append(d)
        found = 1
        break
    l = os.path.join(l , "..")

if found == 0:
    print("Cannot find locally compiled _openmeeg.so file")
    print("Please goto openmeeg repository or a subdir of it")
    exit(-1)
