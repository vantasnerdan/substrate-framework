"""NEW tiny scalar/affine-baseline oracle, not an A26 dynamics replay.

Block order: a separately written polynomial -> exact derivatives and offset
identity -> finite full-density-matrix drift -> production mean/probe/energy.
Support is an independent conservative scalar; ionic quadratic curvature is
already in fixed C. All parameters here are mathematical literals, not Ni/H.
"""
from dataclasses import replace

import numpy as np
import sympy as sp


def polynomial():
    """Encode scalar entries directly, without calling any production action."""
    a,b,c,d,u,v,w,p,q,t = sp.symbols('a b c d u v w p q t', real=True)
    Q = sp.Rational
    # rho=[[a,b+ic],[b-ic,d]], two grid observables, three ionic coordinates.
    n = sp.Matrix([Q(3,5)*a+Q(1,5)*b+Q(3,10)*d,
                   Q(1,5)*a-Q(3,25)*b+Q(1,10)*c+Q(7,10)*d])
    gradient = sp.Matrix([Q(1,5)*a+Q(7,50)*b-Q(3,50)*c-Q(1,10)*d,
                          -Q(2,25)*a+Q(1,25)*b+Q(3,25)*c+Q(9,50)*d,
                          Q(11,100)*a-Q(1,10)*b+Q(1,25)*c+Q(7,100)*d])
    delta = n+sp.Matrix([Q(3,100),Q(7,100)])-sp.Matrix([Q(21,100),Q(31,100)])
    K = sp.Matrix([[Q(4,5),Q(7,50)],[Q(7,50),Q(3,5)]])
    C = sp.Matrix([[Q(6,5),Q(9,100),-Q(1,25)],[Q(9,100),Q(9,10),Q(3,50)],[-Q(1,25),Q(3,50),Q(11,10)]])
    displacement = sp.Matrix([u,v,w])
    ionic = sp.Matrix([Q(19,100),-Q(13,100),Q(1,20)])
    held = sp.Matrix([-Q(3,100),Q(1,20),Q(1,100)])
    support = sp.Matrix([Q(1,50),Q(1,100),-Q(1,100)])
    base = Q(2,5)*a+Q(17,50)*b-Q(6,25)*c-Q(1,5)*d
    electron = base+Q(17,20)*(delta.T*K*delta)[0]
    energy = (electron+(displacement.T*(gradient+held-ionic-support))[0]
              +(displacement.T*C*displacement)[0]/2+p*p/Q(54,10)+q*q/Q(102,10)+t*t/Q(72,10))
    return dict(symbols=(a,b,c,d,u,v,w,p,q,t), energy=energy, electron=electron,
                gradient=gradient, delta=delta, C=C, ionic=ionic, held=held,
                support=support, displacement=displacement)


def check_polynomial(exact, exact_mutation):
    f = polynomial()
    a,b,c,d,u,v,w,p,q,t = f['symbols']
    for axis, coordinate in enumerate((u,v,w)):
        exact('scalar_force_'+str(axis),-sp.diff(f['energy'],coordinate)
              -(f['ionic']-f['gradient']-f['held']-f['C']*f['displacement']+f['support'])[axis])
        exact_mutation('electronic_force_sign_'+str(axis),2*f['gradient'][axis])
    exact('fixed_joint_curvature',sp.hessian(f['energy'],(u,v,w))-f['C'])
    gbar = sp.Matrix([sp.Rational(3,25),-sp.Rational(3,50),sp.Rational(3,100)])
    fbar = f['ionic']-gbar  # DEFINED by same scalar, never independently fitted.
    old_linear = (f['displacement'].T*(f['gradient']+f['held']-gbar-fbar-f['support']))[0]
    new_linear = (f['displacement'].T*(f['gradient']+f['held']-f['ionic']-f['support']))[0]
    exact('common_affine_baseline',old_linear-new_linear)
    exact_mutation('omit_old_gbar_restoration',(f['displacement'].T*gbar)[0])
    exact_mutation('processed_total_plus_electronic_double_count',(f['displacement'].T*f['gradient'])[0])
    extra_ionic = sp.Matrix([[sp.Rational(1,3),sp.Rational(1,20),0],[sp.Rational(1,20),sp.Rational(1,4),0],[0,0,sp.Rational(1,5)]])
    wrong = f['energy']+(f['displacement'].T*extra_ionic*f['displacement'])[0]/2
    exact_mutation('displaced_ewald_quadratic_double_count',sp.hessian(wrong,(u,v,w))-f['C'])
    exact_mutation('unsupported_state_curvature',sp.hessian(f['energy']+a*u*u,(u,v,w))-f['C'])
    exact('support_stays_separate',sp.diff(f['energy']+(f['displacement'].T*f['support'])[0],u)
          -sp.diff(f['energy'],u)-f['support'][0])


def finite_model():
    """Candidate API adapter for declared literal matrices; no comparator RHS."""
    from substrate_framework.physical_factor_preparation import ActionFactors, FactorFunctional, FactorPreparation
    from substrate_framework.prepared_reference import PreparedReference, PreparedSector
    from substrate_framework.prepared_factor_dynamics import PreparedFactorDynamics
    h0 = np.array([[.4,.17-.12j],[.17+.12j,-.2]])
    density = np.array([[[.6,.1],[.1,.3]],[[.2,-.06+.05j],[-.06-.05j,.7]]])
    g = np.array([[[.2,.07-.03j],[.07+.03j,-.1]],
                  [[-.08,.02+.06j],[.02-.06j,.18]],
                  [[.11,-.05+.02j],[-.05-.02j,.07]]])
    K = np.array([[.8,.14],[.14,.6]])
    rho = np.array([[.43,.08+.04j],[.08-.04j,.27]])

    class Coupling:
        npw, nr, ndof = 2, 2, 3
        grid_measure_bohr3, cell_volume_bohr3 = 1.7, 3.4

        def moments(self, left, right, coefficient=None):
            matrix = left@right.conj().T if coefficient is None else left@coefficient@right.conj().T
            return np.einsum('rij,ji->r',density,matrix),np.einsum('aij,ji->a',g,matrix)

        def potential_action(self, field, vectors):
            return 1.7*np.einsum('r,rij->ij',field,density)@vectors

        def ionic_action(self, direction, vectors):
            return np.einsum('a,aij->ij',direction,g)@vectors

    psi = np.eye(2,dtype=complex)
    action = ActionFactors(psi,psi,'I','I','1','composition-original','literal finite identity')
    functional = FactorFunctional.from_actions(psi,{'I':action},factor_identity='composition-raw',pw_basis_identity='composition-original')
    preparation = FactorPreparation(functional,rho,'explicit independent finite quasifree mode')
    sector = PreparedSector(preparation,psi,Coupling(),lambda x:h0@x,0,1.,'composition-original','literal Hermitian matrices')
    reference = PreparedReference(
        nbar=np.array([[.21],[.31]]),ionic_force=np.array([.19,-.13,.05]),
        held_density=np.array([[.03],[.07]]),held_gradient=np.array([-.03,.05,.01]),
        curvature=np.array([[1.2,.09,-.04],[.09,.9,.06],[-.04,.06,1.1]]),
        masses=np.array([2.7,5.1,3.6]),support_force=np.array([.02,.01,-.01]),
        kernel=lambda n:K@n,ncells=1,expansion_authority='fixed joint-second-order polynomial, no extra ionic curvature',
        ionic_authority='explicit finite linear ionic scalar, not material Ewald',
        coverage_authority='one complete two-mode space plus explicit held density/gradient',
        support_authority='independent constant conservative support, not cancellation',
        coordinate_identity='three literal Cartesian coordinates',state_authority='finite nonstationary quasifree density')
    model = PreparedFactorDynamics(reference,(sector,))
    state = model.initial_mean(np.array([.11,-.07,.04]),np.array([.08,-.13,.06]))
    return model,state,rho


def check_composition(check,mutation):
    from substrate_framework.prepared_covariance_source import total_force_probe,probe_mean
    from substrate_framework.prepared_factor_dynamics import PreparedFactorDynamics
    from substrate_framework.prepared_time_observables import boundary_diagnostics
    f = polynomial()
    variables = f['symbols']
    values = (.43,.08,.04,.27,.11,-.07,.04,.08,-.13,.06)
    substitutions = dict(zip(variables,values))

    def numeric(expression):
        return np.asarray(sp.matrix2numpy(expression.subs(substitutions),dtype=complex)) if isinstance(expression,sp.MatrixBase) else complex(expression.subs(substitutions))

    # Differentiate the comparator polynomial, not the candidate force/RHS.
    force = numeric(sp.Matrix([-sp.diff(f['energy'],x) for x in variables[4:7]])).ravel()
    velocity = numeric(sp.Matrix([sp.diff(f['energy'],x) for x in variables[7:10]])).ravel()
    a,b,c,d = variables[:4]
    off = (sp.diff(f['energy'],b)+sp.I*sp.diff(f['energy'],c))/2
    h = numeric(sp.Matrix([[sp.diff(f['energy'],a),off],[sp.conjugate(off),sp.diff(f['energy'],d)]]))
    model,state,rho = finite_model()
    evaluation = model.evaluate_mean(state)
    drift = -1j*(h@rho-rho@h)
    psi,dot = state.frames[0],evaluation.rhs.frames[0]
    candidate_drift = dot@rho@psi.conj().T+psi@rho@dot.conj().T
    check('polynomial_energy',evaluation.energy,numeric(f['energy']))
    check('polynomial_electronic_variation',evaluation.h_frames[0],h)
    check('nonstationary_electronic_drift',candidate_drift,drift)
    mutation('falsely_stationary_electrons',np.zeros_like(drift),drift)
    check('polynomial_force',evaluation.rhs.p,force)
    check('polynomial_acceleration',evaluation.acceleration,force/model.reference.masses)
    mutation('fitted_zero_acceleration',np.zeros(3),force/model.reference.masses)
    check('polynomial_velocity',evaluation.rhs.u,velocity)
    delta = numeric(f['delta']).reshape(2,1)
    check('explicit_held_density_mismatch',evaluation.delta_density,delta)
    mutation('assumed_reference_density',np.zeros_like(delta),delta)
    mutation('omitted_held_density',evaluation.delta_density-model.reference.held_density,delta)
    mutation('omitted_held_gradient',force+model.reference.held_gradient,force)
    check('boundary_energy',boundary_diagnostics(model,state)['energy'],numeric(f['energy']))
    direction = np.array([.7,-.3,.2])
    probe = total_force_probe(model,direction,identity='full affine scalar force')
    check('affine_probe_mean',probe_mean(model,state,probe),direction@force)
    mutation('missing_affine_constant',probe_mean(model,state,replace(probe,constant=0)),direction@force)
    without_support = PreparedFactorDynamics(replace(model.reference,support_force=np.zeros(3),support_authority='explicit absent support'),model.sectors)
    no_support = without_support.evaluate_mean(state)
    check('separate_support_force',evaluation.rhs.p-no_support.rhs.p,model.reference.support_force)
    check('separate_support_energy',evaluation.energy-no_support.energy,-state.u@model.reference.support_force)
    mutation('cancellation_support',force-force,force)
    mutation('processed_total_baseline_double_count',force-evaluation.active_gradient,force)
    extra = np.diag([1/3,1/4,1/5])
    mutation('extra_ionic_curvature_force',force-extra@state.u,force)
