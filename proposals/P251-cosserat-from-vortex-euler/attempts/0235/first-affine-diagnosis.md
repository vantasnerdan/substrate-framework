# First affine oracle diagnosis

The first execution passed the three X-axis predicates and then exposed
an incorrect extension of the raw-lift identity to Y/Z. Source SHA256
was `367872716b6ec91c4fcb8f42688828e21456568cae1aedff2c16ce51be81c5cd`;
the complete failing stdout is preserved as first-affine.stdout, hash
`274a9b75de8afa858e66e0c0de541f8193b2f99847a4dcb7f0258783005ca83b`.

For ALL axes the raw q_e=A_e u-partial_e u=-grad u_e. The actual C015
negative-helicity return is +partial_e u, so the COMPLETE first velocity
is A_e u and the complete material-rate row is zero. The return vanishes
only on X. C016's z forcing is -partial_e u; combined with q_e, its full
physical affine velocity is still R_e+t(w_R-2partial_e u), exactly as
obtained by independent Cartesian Euler substitution. Thus the physical
three-axis/current calculation survives, while the source bookkeeping
and its falsely simplified first predicates are corrected explicitly.

The repair changes those raw-versus-complete rows in proof and verifier.
It also makes the pressure mutation test compute the actual omitted-term
Euler residual instead of merely checking a nonzero source derivative.
This is the same bounded oracle's correction, not a new successful run
of the old incorrect predicates. No whole-run success is assigned to the
preserved first execution.
