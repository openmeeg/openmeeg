# Domain Description 1.1

MeshFile "HeadNNc4.vtp"

Interfaces 5 Interface

Interface North: +north -cut 
Interface South: +south +cut
Interface Cortex: +north +south
Interface Skull: skull
Interface Scalp: 5

Domains 5

Domain NORTH: -North
Domain SOUTH: -South
Domain SKULL: -Skull +Cortex
Domain SCALP: -Scalp +Skull
Domain AIR: +Scalp
