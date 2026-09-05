# Actual cell responses, not magnetic elastic coefficients

The parent registered this expansion and validated its proposal (263/12)
before source-body inspection. All calculations use physical Euler time.

For a steady cell u with pressure p, let t=-partial_U u and q=-partial_U p.
The actual first slow corrector for n perpendicular to U solves

    (u.grad)W+(W.grad)u+grad Pi=-(u.n)t-n q,
    div W=-n.t,  <W>=0.
    c^2 U=-P_n <u tensor W+W tensor u>n.

These are the full pressure and mean equations from 0129, not an energy
restricted to a selected displacement. A particular solution alone does
not establish uniqueness modulo the full stationary Euler kernel.

## Separable planar cellular candidate

Take u=(-b sin y,a sin x,0), p=-ab cos x cos y, n=e_x, U=e_y.
The exact corrector is W=(0,-b sin y,0), Pi=ab sin x sin y.
Direct differentiation gives the stated divergence and full momentum
equations; the response is c^2 U=-(b^2/2)e_y. Interchanging axes gives
the coefficient -a^2/2. This is an exact cell-response calculation, not
yet a spectral theorem for the cellular flow's unrestricted kernel.

The apparent alternative using the first Arnold criterion Delta psi=f(psi),
f'>0 does not yield a nonconstant smooth periodic cell: differentiating
and integrating yields

    -integral |grad partial_i psi|^2
      =integral f'(psi)|partial_i psi|^2.

It forces each derivative to vanish. This refutes that **periodic first-
criterion route**, not Arnold stability with boundaries, a different
Casimir criterion, or every steady Euler cell.

## Exact three-dimensional first-shell calculation

The candidate u=(sin z+C cos y,sin x+cos z,C sin y+cos x) has curl u=u.
For C=1 and C=0, n=e_z, U=e_x, the exact corrector W=(cos x,0,0)
solves the cell equation with its retained gradient pressure. The response
is -e_x/2. The script solves all 18 real first-shell corrector coefficients;
all free coefficients cancel from this response. It evaluates every
Fourier product, rather than projecting the residual onto that shell.
For diagonal n=(1,1,0) and (1,1,1), U=(1,-1,0), the same finite ansatz is
inconsistent by exact augmented rank. That is a representation-scoped
verdict; a broader corrector space remains available.

## Primary-source boundary and next candidate

Moffatt, JFM166 (1986), 359--378, primary
http://www.damtp.cam.ac.uk/user/hkm2/PDFs/Moffatt_1986_JFM_166_359.pdf,
is archived by the parent at /tmp/moffatt-pr199-1986.pdf (SHA256
c518d49675c128d7327666535cc3eed5d3affa390df512456b89c8c9d704ffc7).
Its section 7 explicitly says its proposed mean dynamical equation (7.26)
is not rigorously established. The omitted small-velocity/large-vorticity
coupling is exactly why the full pressure corrector above is retained.
Its magnetic elastic waves are not imported as Euler waves.

The failure-derived array candidate now uses actual opposite-circulation
streets. Tkachenko, JETP37 (1973), 1148--1149,
https://www.jetp.ras.ru/cgi-bin/dn/e_037_06_1148.pdf, supplies a primary
point-vortex oscillation comparator, not smooth-core or positive-energy
authority. Garcia, Nonlinearity33 (2020), 1625--1676,
https://arxiv.org/abs/1905.07728, constructs traveling vortex *patches*;
its abstract alone is not used as a smooth-vorticity or stability theorem.
The street calculation below is independently derived from Biot--Savart.
