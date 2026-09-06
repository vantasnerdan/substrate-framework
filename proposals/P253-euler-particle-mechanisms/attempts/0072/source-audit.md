# Source and authority audit

This attempt is an exact analytic continuation of the independently reviewed
0069 escape route. It uses no new empirical comparator and imports no theorem
claim from an unreviewed paper. The homogeneous Fourier transform is derived
as a tempered-distribution identity with the stated convention; the Euler
tail evolution follows directly from the Newton pressure representation and
the declared weighted classical domain.

The reviewed 0069 result supplies only the existence and checked algebra of
the representative (12) and its leading cross kernel. This attempt rederives
those formulas inside the complete harmonic classification. Accepted Euler
local well-posedness is used only to name the classical time interval; the
weighted tail preservation is proved by the displayed pressure split rather
than attributed to a source.

No compact-vorticity moment theorem is imported: the whole point of this
route is that `curl u` is not `L^1`. No KKS, charge, detector, quantum, or
relativistic conclusion is inherited from earlier attempts.

## Primary stationary-homogeneous source

Roman Shvydkoy, *Homogeneous solutions to the 3D Euler system*,
arXiv:1510.03378v1, was fetched from
`https://arxiv.org/pdf/1510.03378v1` after activation. Cached PDF:

    /tmp/primary-source-cache/P253-0072/Shvydkoy-1510.03378v1.pdf
    SHA-256 71b277e2b78c2e8a1d14994c7267458795ee630bd555e0cb856fcc51402f5706

Page 3, equations (3a)--(3c), gives for
`V=r^(-alpha)(v+f n)`

    (2-alpha)f+div_S v=0,
    v dot grad_S f=|v|^2+alpha f^2+2 alpha p,
    (1-alpha)f v+nabla_v v=-grad_S p.

At `alpha=2` these are exactly (3) and (17a). Proposition 5.2 on page 12
classifies only the axisymmetric `alpha=2` case: the radial constant-flux
solution is the sole member, and it is excluded by the smooth source-free
flux row here. The paper explicitly says the general case remains open. It
therefore supports the sphere equations and the axisymmetric boundary, not a
general nonaxisymmetric degree-minus-two no-go. The fixed-frame texture is
nonaxisymmetric and is tested directly only at the stated stress and
`l=1`-projection scopes.
