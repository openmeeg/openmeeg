file(GLOB files "Head*")
MESSAGE("REMOVE ${files}")  # Just verbose in case is needed
file(REMOVE ${files})