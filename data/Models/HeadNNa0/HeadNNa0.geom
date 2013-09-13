# Domain Description 1.1

MeshFile "HeadNNa0.vtp"

Interfaces 3

Interface Cortex: north south
Interface Skull: skull
Interface Scalp: scalp

Domains 4

Domain SOUTH: -Cortex
Domain SKULL: -Skull +Cortex
Domain SCALP: -Scalp +Skull
Domain AIR: +Scalp
