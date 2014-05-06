NEURON {
POINT_PROCESS Gap
RANGE vgap
RANGE g, i
NONSPECIFIC_CURRENT i
}
PARAMETER { g = 1 (nanosiemens) }
ASSIGNED {
v (millivolt)
vgap (millivolt)
i (nanoamp)
}
BREAKPOINT { i = (v - vgap)*(g*1e-3) }
