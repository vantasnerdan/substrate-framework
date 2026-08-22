verifier: verify_els008
exit_code: 0
verdict: ALL 16 CHECKS PASS [C-ELS-008]

C-ELS-008 (barotropic stiff EOS, p = K*rho, premise P0): the filtered-balance
rung extends to the compressible closure with the equation of state as a named,
independently falsifiable premise. Exact symbolic results: lambda_total =
rho0*K + E_f*L_v/15 and mu_tangle = E_f*L_v/15 by coefficient matching;
Christoffel speeds c_P^2 = K + E_f*L_v/(5 rho0), c_S^2 = E_f*L_v/(15 rho0);
fluid limit mu->0 restores c_P^2 = K with no shear branch; Poisson ratio
interpolates nu = 1/2 (EOS only) to 1/4 (tangle only), making nu an
EOS/tangle discriminator. Dynamics companion: longitudinal oscillator at
c_P matches the symbolic speed; half-strength bulk modulus oscillates
measurably off-shell. Mutations: dropped body force breaks the solution
class; sign-flipped closure loses strong ellipticity.
