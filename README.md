# Effective_Mass
Calculate effective masses of electron, hole, and exciton from band structure

Calculate *nscf* with Quantum Espresso using 4 k-points:

```
K_POINTS tripa
 4
 0.0    0.0    0.0    1.
 0.02   0.     0.     1.
 0.     0.02   0.     1.
 0.     0.     0.02   1.

```

and select necessary bands for study.
