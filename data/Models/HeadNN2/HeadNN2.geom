# Domain Description 1.1

MeshFile "HeadNN2.vtp"

Interfaces 4 Interface

Interface North: +north +cut 
Interface South: +south -cut
Interface Skull: north south skull
Interface Scalp: skull scalp

Domains 5

Domain NORTH: -North
Domain SOUTH: -South
Domain SKULL: -Skull +North +South
Domain SCALP: -Scalp +Skull
Domain AIR: +Scalp
