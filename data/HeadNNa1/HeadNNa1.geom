# Domain Description 1.1

Meshes 5
Mesh cut: "cut.tri"
Mesh north: "north.tri"
Mesh south: "south.tri"
Mesh skull: "skull.tri"
Mesh scalp: "scalp.tri"

Interfaces 5

Interface North: +north +cut 
Interface South: +south -cut
Interface Cortex: north south
Interface Skull: skull
Interface Scalp: scalp

Domains 5

Domain NORTH: -North
Domain SOUTH: -South
Domain SKULL: -Skull +Cortex
Domain SCALP: -Scalp +Skull
Domain AIR: +Scalp
