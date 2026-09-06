# P253/0015 bounded normalization and compact-extension correction

Date: 2026-09-06
Origin: P253/0023 single bounded review finding

The symbols \(E,P,\Omega,\mathscr A\) in `construction.md` are now declared
normalized. Equations (8), (9), (16), and (18) restore the common physical
factor \(2\pi\rho_m\), where \(\rho_m\) is material mass density and is not
the canonical coordinate \(\rho=r^2/2\):

\[
 E_{\rm kin}^{\rm phys}=2\pi\rho_mE,\quad
 I_z^{\rm phys}=2\pi\rho_mP,\quad
 \Omega^{\rm phys}=2\pi\rho_m\Omega,\quad
 \mathscr A^{\rm phys}=2\pi\rho_m\mathscr A.
\]

Hence
\(i_X\Omega^{\rm phys}=d\mathscr H_c^{\rm phys}\) is equivalent to the
normalized Hamilton equation, and the positive common factor leaves all
criticality, Jacobi, indefiniteness, and route-verdict statements unchanged.

Equation (7a) also chooses a compact cutoff Hamiltonian equal to
\(\Psi-c_\epsilon\rho\) on a neighborhood of the compact union of periodic
patch supports, which lies away from the axis. Its compactly supported
area-preserving flow agrees with the translating-frame transport on both
patches, supplying the claimed `Diff_c` leaf membership without treating the
uniform translation at infinity as compactly supported.

Hash receipt:

- `construction.md` before correction:
  `1d146f6e1d68175a3844e0d27ea3acda6bbe41bb85e88af93905d65fb2a38eb6`
- `construction.md` after correction:
  `90eec03cd419ca92e6796e6e4d3839b6c6d2152ffa2b04ef49048eacd3a5bff9`
- unchanged `action-sign-correction.md`:
  `a64811d54db5d2c8348e1fa6ab2098a8f886f189290c4067014772c3c1546d80`
- unchanged `verify_exact_structure.py`:
  `0e9350ef0b7b0b7c1eca869096395c8b92cd9e429522466eac98c7fe406cf981`
- unchanged `exact-check.stdout`:
  `04c84f52c49e3d10b666bb3c09cf19a4652bda1a69e40bff8dbd3ddb50086287`
- unchanged `exact-check.exit`:
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`

Per the bounded-check instruction, the unchanged exact oracle was not rerun.
