# Positive nonunit spin overlap preserves the standard physical-angle chart

This is a conditional interface for the actual stationary mode sought in
0181. It does not replace its mode/current construction. Let rho,j,gap>0,
gap=omega², and let eta>0 be the measured stationary spin overlap. Keep
the literal acoustic tag inertia I, the independent acoustic gradient
mass m_a and optical gradient mass j2. In a transverse curl helicity h,
h²=k², use the SAME uneliminated two-field reference as0172:

    M0=diag(rho+m_a*k²,j+j2*k²),
    K0=diag(rho*a*k²,j*gap+(j*B+gap*j2)*k²),
    U=(1-I*k²/(4rho))X-eta*j*h*q/(2rho),
    Phi=q+h*X/2.

The actual nonunit overlap creates mixed physical mass
M_U,Phi=(eta-1)*j*h/2. This is not deleted or declared a zero overlap.
Instead define the explicit current-improved displacement

    U_c=U+(eta-1)*j*h*Phi/(2rho), Phi_c=Phi.

Composing the actual maps gives exactly, through the second jet,

    U_c=[1-I_eff*k²/(4rho)]X-j*h*q/(2rho),
    I_eff=I+(1-eta)*j.

Thus it is the original complete0172 map with its acoustic gradient mass
changed, not with its physical angle rescaled. In particular a uniform
rigid rotation has constant Phi, curl Phi=0: the correction leaves U and
the unit angular response unchanged. I_eff is an algebraic coefficient,
not a replacement for the literal positive tagged moment I; I is still
in the physical observation.

Pull back BOTH M0 and K0 through this composed map. The physical gradient
masses become

    m_U=m_a+I_eff/2-j/4,
    m_Phi=j2-j²/(4rho).

The usual derivative map diag(1-m_U*k²/(2rho),
1-m_Phi*k²/(2j)) normalizes both masses while transforming the potential
at the same time. The resulting standard coefficients are still

    mu=rho*a, alpha=j*gap/4,
    C_T=j*B-j²*gap/(4rho), C_L=j*B_L.

Their sign condition is the same finite positive density inequality as
0172. No eta-dependent fitted modulus appears. The actual measured
displacement is reconstructed as

    U=U_c-(eta-1)*j*curl(Phi)/(2rho),

followed by the displayed gradient normalization. For the optical branch,
X=0 in the diagonal reference and Phi=q, so the LITERAL physical transfer
is

    U/Phi=-eta*j*h/(2rho)+O(k³).

It is nonzero for every eta>0. At eta=0 it vanishes even though the
current-improved canonical equations retain their usual coupling: that
limit exposes why a coordinate coupling is not by itself a physical
transfer. Here the actual measured row supplies the needed distinction.

All statements above assume the stationary mode's actual j,eta and
integrated spin row. For time-dependent values c(t)=(eta(t)-1)j(t)/(2rho),
the map has U_c,t=U_t+c curl(Phi_t)+c_t curl(Phi), and the pulled-back
action retains c_t and c_tt terms. A painted marker does not acquire a
stationary current by applying this normalization. Its construction
remains0181's task, now with every strictly positive overlap admissible.

Route verdict: established as an exact same-action second-jet interface.
Evidence scope: conditional physical-current normalization, not existence
of the stationary Euler mode or completion of the coupled campaign.
