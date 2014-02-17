# Domain Description 1.1

MeshFile "HeadNNb2.vtp"

Interfaces 5

Interface SphereNorth: spherenorth
Interface SphereSouth: spheresouth
Interface Cortex: cortex
Interface Skull: skull
Interface Scalp: scalp

Domains 6

Domain SPHERENORTH: -SphereNorth
Domain SPHERESOUTH: -SphereSouth
Domain BRAIN: -Cortex +SphereNorth +SphereSouth
Domain SKULL: -Skull +Cortex
Domain SCALP: -Scalp +Skull
Domain AIR: +Scalp
