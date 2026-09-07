# Postactivation momentum-domain correction receipt

## Chronology and byte boundary

P253/0107 was centrally registered and schema-activated with exit `0` on
README SHA-256
`dc8aa8efb88f214b897dbf95bf42a48551fb54e0e3d521be6884e68e07630e0c`.
The activation receipts are preserved byte-for-byte:

- command: `5469bc1a4e263e5211cc04b02d369ba48f98d4da7a5d5cd5a40c24476426451f`;
- stdout: `d7a6df70d40116eccbc5ef95222bcd39f56874086284528ed1e49c9c2cd6676f`;
- stderr: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`;
- exit: `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

Before any Route-A witness conclusion, a coupled-balance audit exposed that
the activated README's velocity-integral momentum was not an invariantly typed
row.  The corrected README SHA-256 is
`59cfe8d5ec379fdded8f34ecdf951ae9a66c6b910603f8dbef1fc1e94f91f345`.
Body conclusions remain paused until central schema replay on that hash.

## Exact calculation that forced the correction

Use

    q=omega-omega_g,
    I_h(omega)=(rho_m/2) integral x cross omega dx,
    f=rho_q E+j cross B.

For a localized forced incompressible Euler solution,

    partial_t omega=curl(u cross omega)+(1/rho_m)curl f.

The component identity, justified first on balls and then by the stated decay,
is

    (1/2) integral x cross curl A dx=integral A dx.

Moreover

    u cross omega=grad(|u|^2/2)-(u dot grad)u,

whose integral vanishes for a localized divergence-free velocity.  Therefore

    d I_h/dt=integral f dx.                               (R1)

With Maxwell equations

    div(epsilon_EM E)=rho_q,
    div B=0,
    epsilon_EM partial_t E=(1/mu_EM)curl B-j,
    partial_t B=-curl E,

put

    sigma_EM=epsilon_EM(E tensor E-|E|^2 I/2)
             +(1/mu_EM)(B tensor B-|B|^2 I/2).

Direct vector calculus gives

    partial_t(epsilon_EM E cross B)
       =div sigma_EM-f.                                  (R2)

Thus, with momentum flux `-sigma_EM`, the expanding-ball limit and vanishing
relative flux give

    d P_EM/dt=-integral f dx,
    P_EM=epsilon_EM integral E cross B dx.               (R3)

Subtracting the traveling carrier balance yields the conserved finite relative
translation moment map

    J_ren=[I_h(omega)-I_h(omega_g)]
          +epsilon_EM integral(E cross B-E_g cross B_g)dx,
    d J_ren/dt=0.                                        (R4)

By contrast, if `v=u-u_g` is divergence-free and belongs to `L1`, Fourier
continuity gives `integral v=0`.  That fact makes the velocity integral a
domain diagnostic, not a second momentum.  Lorentz exchange can change
`I_h`, generate an `O(|x|^-3)` Hodge tail, and take `v` out of `L1` unless a
separate zero-net-force theorem is proved.

## Contract effect

The corrected contract therefore:

1. replaces velocity-integral momentum by `J_ren` as the invariant row;
2. retains `Delta I_h=0` only as an initial carrier/modulation slice, not as an
   independently conserved coupled quantity;
3. stages the row solve: impulse/center first, then `J_ren`;
4. uses the exact lower-triangular derivative

       [[H,0,0],
        [0,C,0],
        [K_I,K_C,G]],

   because displacement columns also change the advected tag and its
   Gauss-reconstructed longitudinal electric field; and
5. preserves the packet-projection order and every frozen Route B/C target.

No H/C/G invertibility, exact fixed-leaf curve, nonlinear contradiction, P2,
or parent conclusion is claimed by this receipt.  Those are the achievements
unlocked after corrected schema replay.
