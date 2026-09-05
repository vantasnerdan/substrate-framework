# 0266 source and exact-oracle receipt

Original pre-centering cellwise-phase.md SHA256 (checkpoint3ae5f69):
6cf3dea30789a39c02f910fcb60ba01ab7fc85e911b5fc62e6d22ce332308794.

Original uncentered verify_material_tensor.py SHA256:
180e49da24de7f2ecc238a4778e8a6fba50c809079cfec59c75d531ac772779a.

The active Cartesian oracle computes all18 symmetric tensor entries by
integrating the literal material expression, derives its oblique transverse
projection and sphere moments, and exposes deletion of the mixed entries.
material-tensor.stdout/exit is exit0. It finds the averaged coefficient
4 i pi^2 C V^2(R^2+s^2)/(15s), not the earlier axis-only scalar.

The first execution stopped on a structural SymPy equality between
different factorizations of the same polynomial. Original script/output
and exit1 are preserved with the _initial/-initial suffixes. Replacing
structural equality by simplifying the exact difference resolves the
implementation failure; no equation, expected coefficient or tolerance
was changed. Targeted Ruff passes. The script ran in Herdr shell w3:p2.

0267 identified the missing central terms in the all-18-entry description.
The active source and oracle now retain both
-deltaX_j A_0,il-deltaX_l A_0,ij, with the actual integrated baseline spin
and centroid variation independent. The baseline shell first moment is
derived separately. This changes the chiral B_yzx entry by
-a_spin delta_X_z and leaves the even coefficients and sphere gain
unchanged. The prior successful uncentered script/output are preserved
with the _uncentered/-uncentered suffixes. The active check and0267's
independent rerun both pass exit0; its bounded correction review accepts
the repair as stated.

Current cellwise-phase.md SHA256:
e5c3e641c9f9614875e57861981182bdc80d3f4bb0d5e1b8b617115a8ff88765.

Current verify_material_tensor.py SHA256:
22002752fcee1e333f2a23b4966bf45058e8ec22b725cbd37ea384fc46b37e0e.

The oracle corroborates the explicit circular leading tensor. It does
not alone prove cellwise-source preparation, finite-window pressure
transfer, centroid-law cancellation or finite-R persistence; those are
the source-level0267 independent review boundary. No parent completion
or accepted claim is inferred from the tally.
