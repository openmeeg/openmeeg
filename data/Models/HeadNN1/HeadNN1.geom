# Domain Description 1.1

MeshFile "HeadNN1.vtp"

Interfaces 3 Interface

Interface Brain: Cortex
Interface Skull: Skull
Interface Scalp: Scalp

Domains 4

Domain Scalp: Skull -Scalp
Domain Brain: -Brain
Domain Air: Scalp
Domain Skull: Brain -Skull

