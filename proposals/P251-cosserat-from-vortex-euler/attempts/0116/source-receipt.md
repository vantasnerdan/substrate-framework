# Primary-source and scale receipt

Read before source-dependent selection: Enciso, Peralta-Salas and Torres
de Lizaur, [Knotted structures in high-energy Beltrami fields on the torus
and the sphere](https://arxiv.org/pdf/1505.01605), Theorems 1.1, 2.1, 2.2
and the torus construction in Section 7; [published version](https://numdam.org/articles/10.24033/asens.2337/).

The torus theorem supplies smooth exact curl eigenfields with prescribed
knotted tubes in a contractible region. Its inverse-localization statement
approximates a Euclidean eigenfield in any fixed finite differentiability
norm after rescaling, for sufficiently large odd eigenvalue. The torus
construction uses finitely many integer Fourier modes on one curl shell.
Its approximation and structural-stability statements concern the local
field and tube geometry; they do not assert an optical Bloch band or a
coarse constitutive law.

Archived source paths are `../../sources/1505.01605.pdf` and
`../../sources/1505.01605.pdf.txt`.
SHA256, respectively:

```
8142ed18ae0d7d11164abed4115405be87f00815b585192ce44aaf87218a4c6c
a6e4e12b3a29d4a0cbe6e1663f7d759d67543a15d40abd4d22209bab7ba76bf7
```

The following scaling is derived here, not an additional source theorem.
On `(R/2πZ)^3`, if `curl u_Λ=Λu_Λ` and `u_Λ(y/Λ)` approximates the
chosen unit-curl local field, core lengths and traversal times at bounded
velocity scale as `Λ^-1`, while the local angular rate scales as `Λ`.
For a target eigenvalue of either sign, first rescale the Euclidean model
to unit absolute eigenvalue and choose the matching helicity. Fix the
periodic background and its finite norm bounds before choosing packet
carrier, Bloch wave number, or nonlinear perturbation amplitude.

The local packet recursion is imported from 0112/0114 with its fixed
physical-frame scope. Transferring it uses the periodic Leray operator
directly, as explained in `periodic-bloch-transfer.md`; local field
approximation alone is not a transfer theorem for a nonlocal projector.

Classical smooth Euler local evolution is used only on a fixed finite
supercell and time interval, with perturbation amplitude chosen last.
For provenance, the primary compact-manifold framework is Ebin and
Marsden, [Groups of diffeomorphisms and the motion of an incompressible
fluid](https://annals.math.princeton.edu/1970/92-1/p04). No global-in-time
regularity or uniform-in-supercell lifespan is claimed; the perturbative
finite-time estimate needed here is displayed in the proof.
