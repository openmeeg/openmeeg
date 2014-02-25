# Domain Description 1.1

Interfaces 3

Interface Skull: "skull.vtk"
Interface Cortex: "brain.vtk"
Interface Head: "head.vtk"

Domains 4

Domain Scalp: Skull -Head
Domain Brain: -Cortex
Domain Air: Head
Domain Skull: Cortex -Skull
