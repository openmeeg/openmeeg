# Domain Description 1.1

Interfaces 4

Interface scalp: "scalp.off"
Interface skull: "skull.off"
Interface left:  "left.off"
Interface right: "right.off"

Domains 5

Domain Air: scalp
Domain Scalp: skull -scalp
Domain Left: -left
Domain Right: -right
Domain Skull: +left +right -skull
