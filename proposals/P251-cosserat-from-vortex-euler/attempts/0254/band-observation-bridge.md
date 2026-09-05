# A fixed open frequency band and the actual observation error

## 1. Constructive finite-window density

Fix T>0, m>=0, an open real interval J and omega0 in J. Every complex
f in C^m([-T,T]) can be approximated in that norm by a finite linear
combination of exp(i omega t) with omega in J. To construct it, approximate
exp(-i omega0 t)f(t) by a polynomial P(t)=sum_(n=0)^N p_n t^n.
Polynomial density in C^m follows by approximating the mth derivative
uniformly and integrating m times, retaining the lower initial jets.

For a nonzero real h sufficiently small that omega0+[0,N]h lies in J,
put g_h(t)=(exp(i h t)-1)/(i h). Then

    F_h(t)=exp(i omega0 t) P(g_h(t))
          =sum_(j=0)^N c_j exp(i(omega0+jh)t),
    c_j=sum_(n=j)^N p_n binom(n,j)(-1)^(n-j)/(i h)^n.

This is a finite frequency construction, with

    sum_j |c_j| <= sum_n |p_n|(2/|h|)^n.                (1)

Since g_h(t)=integral_0^t exp(i h s)ds, the zeroth error is at most
|h|T^2/2. The first derivative error is at most |h|T, and for j>=2,
g_h^(j)=(i h)^(j-1)exp(i h t). Thus g_h -> t in every fixed C^m,
with O(|h|) for |h|<=1. Finite products and the fixed carrier preserve
this bound, so F_h -> exp(i omega0 t)P(t) in C^m. Polynomial error
first, h second, proves the theorem with entirely finite coefficients.

For real histories, take real parts, which uses sine/cosine quadratures
at the same real frequencies. Even/odd projection on[-T,T] gives the
corresponding parity targets. The required physical quadratures or paired
source law must actually exist. No frequency tending to zero is required.
The coefficient norm can be very large; (1) is retained rather than treated
as an innocuous normalization.

The dual proof gives the same scope: a C^m functional annihilating all
these exponentials defines an entire function of complex frequency. Its
vanishing on J implies it vanishes identically, so it annihilates every
polynomial and hence every C^m function. This proof does not furnish an
additional Euler supplier.

## 2. Quasimodes normalized by their physical gain

Let S(t) be the actual strongly continuous group for the full linear state,
with generator L and sup_(|t|<=T)||S(t)||<=M_T. Let C be one fixed bounded
physical scalar observation. For a real omega and u in Dom L^(m+1), set
r=(L-i omega)u and gamma=C u !=0. Duhamel and the difference-of-powers
identity give, for j>=1,

    S(t)u-exp(i omega t)u
      =integral_0^t exp(i omega(t-s))S(s)r ds,
    [L^j-(i omega)^j]u
      =sum_(a=0)^(j-1)(i omega)^(j-1-a)L^a r.

Consequently the jth time derivative of C S(t)(u/gamma)-exp(i omega t)
has sup norm at most

    ||C|| M_T / |gamma| *
      [sum_(a=0)^(j-1)|omega|^(j-1-a)||L^a r||
                                  +T|omega|^j||r||],   (2)

while for j=0 it is at most ||C||M_T T||r||/|gamma|.
A small raw residual is insufficient: its required graph norms must be
small relative to the measured gain. The gain and residual refer to the
same physical operator, observation and state.

If each omega in J has actual smooth states with all ratios in(2) tending
to zero, choose the finite frequencies/weights in section1 first. Then
choose their finite states so each normalized output error is smaller
than the desired remaining accuracy divided by 1+sum_j|c_j|. Their sum
is one actual initial state on the same background. Its finite norm is
bounded by sum_j |c_j| ||u_j||/|gamma_j| and must enter every later
pressure, action, current and long-wave remainder estimate before K is
chosen. Neither the background nor the observation changes during this
linear combination.

## 3. Several independent channels and physical constraints

For q independent outputs a sufficient replacement is a q-column actual
source family at each frequency whose measured gain matrix Gamma is
invertible. Normalize the whole family by Gamma^(-1), retain its norm,
and apply the columnwise version of(2). This gives the same approximation
for vector histories. A single nonzero scalar gain does not prove this
matrix assertion. Exact relations between physical current and displacement
can instead define a smaller compatible target space, but that range must
be derived from the actual observation operator.

An advected/time-dependent detector is not automatically this C. It needs
an actual stationary augmented representation with the stated group bound,
or a separate calculation of its time-dependent response kernels. Changing
a material tag for each packet likewise needs one common positive tag/law
and its full gain matrix before summation; packet-specific observers cannot
be combined as if they were the same detector.

This is the precise0250 achievement unlocked: a real open response band
with controlled normalized generator residuals and the actual required
observation range on ONE fixed background. Passive particle frequencies,
finite Taylor ranks, an unobserved Weyl sequence, or an energy-only phase
normalizer do not satisfy those hypotheses. The parent needs the same
source family and every actual action/current cross form afterwards.
