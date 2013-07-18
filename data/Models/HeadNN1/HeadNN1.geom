# Domain Description 1.1

MeshFile "HeadNN1.vtp"

Interfaces 3 NamedInterface

Brain: Cortex
Skull: Skull
Scalp: Scalp

Domains 4

Domain Scalp Skull -Scalp
Domain Brain -Brain
Domain Air Scalp
Domain Skull Brain -Skull

