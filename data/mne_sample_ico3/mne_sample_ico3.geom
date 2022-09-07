# Domain Description 1.1

Interfaces 3

Interface Skull: "outer_skull.tri"
Interface Cortex: "inner_skull.tri"
Interface Head: "outer_skin.tri"

Domains 4

Domain Scalp: Skull -Head
Domain Brain: -Cortex
Domain Air: Head
Domain Skull: Cortex -Skull
