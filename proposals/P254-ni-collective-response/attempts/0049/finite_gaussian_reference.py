"""Finite product preparations and independent exact full quadratic evolution.

Only oracle producers: not part of the general supplied-kernel solver API.
Boson inputs are physical real Cartesian x, p in hbar=1; no oscillator scaling.
"""
import numpy as np
from scipy.linalg import expm


def fermion_case(times):
    rotation = expm(-1j*np.array([[.2, .7-.3j], [.7+.3j, -.1]]))
    rho = rotation @ np.diag([0., 1.]) @ rotation.conj().T
    bath_rho = np.array([[.31, .08+.03j], [.08-.03j, .64]])
    hs = np.array([[.4, .25+.13j], [.25-.13j, -.3]])
    hb = np.array([[-.2, -.11+.07j], [-.11-.07j, .6]])
    v = np.array([[.37+.08j, .12], [-.09j, .28-.04j]])
    full_h = np.block([[hs,v], [v.conj().T,hb]])
    full_rho = np.block([[rho,np.zeros((2,2))], [np.zeros((2,2)),bath_rho]])
    def components(h, state):
        u = np.array([expm(-1j*h*t) for t in times])
        less = np.array([[1j*a@state@b.conj().T for b in u] for a in u])
        ret = np.zeros_like(less)
        for n,a in enumerate(u):
            for m in range(n+1):
                ret[n,m] = -1j*a@u[m].conj().T
        return ret, less
    source = components(hs,rho)
    br,bl = components(hb,bath_rho)
    kernel = tuple(np.einsum('ab,tsbc,dc->tsad',v,a,v.conj()) for a in (br,bl))
    exact = tuple(a[:,:,:2,:2] for a in components(full_h,full_rho))
    discarded = tuple(a[:,:,:2,:2] for a in components(full_h,np.diag(np.diag(full_rho))))
    return {'source': source, 'kernel': kernel, 'exact': exact,
            'discarded_preparation': discarded, 'rho': rho,
            'rho_eigenvalues': np.linalg.eigvalsh(rho)}


def phase_generator(masses, stiffness):
    d = len(masses)
    return np.block([[np.zeros((d,d)),np.diag(1/np.asarray(masses))],
                     [-stiffness,np.zeros((d,d))]])


def prepared_covariance(masses, frequencies, squeeze):
    masses, frequencies = np.asarray(masses), np.asarray(frequencies)
    d = len(masses)
    # Independent vacuum preparation transformed by a real symplectic shear.
    initial = np.diag(np.r_[1/(2*masses*frequencies), masses*frequencies/2])
    transform = np.block([[np.eye(d),np.zeros((d,d))], [squeeze,np.eye(d)]])
    return transform@initial@transform.T


def boson_components(times, generator, covariance, positions):
    d = len(covariance)//2
    j = np.block([[np.zeros((d,d)),np.eye(d)], [-np.eye(d),np.zeros((d,d))]])
    s = np.array([expm(generator*t)[positions] for t in times])
    less = -1j*np.einsum('tar,rs,ubs->tuab',s,covariance-.5j*j,s)
    ret = np.zeros_like(less)
    for n in range(len(times)):
        # Canonical POSITION coordinates commute at equal times exactly.
        # This is the finite producer's definition, not a solver-side repair.
        for m in range(n):
            ret[n,m] = s[n]@j@s[m].T
    return ret, less


def boson_case(times):
    masses = np.array([1.3,2.1,.7])
    stiffness = np.array([[1.7,.31,.24], [.31,2.8,-.19], [.24,-.19,1.2]])
    cs = prepared_covariance(masses[:2],[.9,1.4], np.array([[.21,.16], [.16,-.14]]))
    # Noncommuting x covariance from a canonical Cartesian rotation.
    angle = .37
    rotation = np.array([[np.cos(angle),-np.sin(angle)], [np.sin(angle),np.cos(angle)]])
    transform = np.block([[rotation,np.zeros((2,2))], [np.zeros((2,2)),rotation]])
    cs = transform@cs@transform.T
    cb = prepared_covariance(masses[2:],[1.1],np.array([[.23]]))
    cov = np.zeros((6,6))
    sys_indices, bath_indices = [0,1,3,4], [2,5]
    cov[np.ix_(sys_indices,sys_indices)] = cs
    cov[np.ix_(bath_indices,bath_indices)] = cb
    source = boson_components(times,phase_generator(masses[:2],stiffness[:2,:2]),cs,[0,1])
    br,bl = boson_components(times,phase_generator(masses[2:],stiffness[2:,2:]),cb,[0])
    coupling = stiffness[:2,2:]
    kernel = tuple(np.einsum('ab,tsbc,dc->tsad',coupling,a,coupling) for a in (br,bl))
    generator = phase_generator(masses,stiffness)
    exact = boson_components(times,generator,cov,[0,1])
    no_cross = np.diag(np.diag(cov))
    discarded = boson_components(times,generator,no_cross,[0,1])
    wrong_mass = boson_components(times,phase_generator(np.ones(3),stiffness),cov,[0,1])
    j = np.block([[np.zeros((3,3)),np.eye(3)], [-np.eye(3),np.zeros((3,3))]])
    return {'source':source,'kernel':kernel,'exact':exact,
            'discarded_preparation':discarded,'wrong_mass':wrong_mass,
            'covariance':cov,'masses':masses,'stiffness':stiffness,
            'uncertainty_min_eigenvalue':float(np.linalg.eigvalsh(cov+.5j*j).min()),
            'coordinate_covariance_stiffness_commutator':float(np.linalg.norm(cs[:2,:2]@stiffness[:2,:2]-stiffness[:2,:2]@cs[:2,:2]))}
