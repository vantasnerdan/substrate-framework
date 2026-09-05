# 0177 — canonical observed Euler response rows

Base v0.177.0/dfe495c. C-CST-013 was independently accepted in0173,
without a correction. Extract its actual first-cell preparation/forcing
and mean acceleration/current rows into additive euler_acoustic functions
calling the existing exact euler_fourier operations.

The canonical API retains the pressure per density, harmonic projection,
transverse macro preparation, two pressure terms and Eulerian/material
current distinction. It returns algebraic rows on supplied actual cell
solutions; it does not solve or certify their PDE evolution. Stationarity,
solenoidality and preparation/domain conditions are checked where exactly
decidable, with symbolic unresolved hypotheses explicit in the interface.

Tests compare independent one-wave Cartesian/covariance expressions,
nonconstant-pressure data exposing both stress and current, exact Euler
stationarity, and invalid domains. Only these new functions and direct
tests are scientific replay targets. The existing acoustic definitions,
Fourier API and all earlier claims remain unchanged. Registry/generated
and memory checks follow the separate promotion transaction. No full
suite is requested or run under the user's blast-radius instruction.
