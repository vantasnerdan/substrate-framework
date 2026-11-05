"""Fixed-original-full-state ordinary-NC electronic onsite third derivative.

Ionic derivatives act ONLY on exp(-ip.tau), never on radial interpolation.
All eight projector product-rule legs survive with their complex phases.
Raw full-column D and every original slot weight belong to the caller. No
spin2, MV refill, occupied-band projection, SCF, electronic response or H/F/J.
Local density includes its full physical 1/Omega through the existing FFT.
"""
from itertools import combinations_with_replacement, permutations
import numpy as np
from scipy.linalg.blas import zgemm
from cubic_position_moments import compressed_overlaps,position_fft


TRIPLES = tuple(combinations_with_replacement(range(3),3))
MOMENTS = tuple(combinations_with_replacement(range(3),degree)
                for degree in range(4))
MOMENTS = tuple(moment for degree in MOMENTS for moment in degree)


def moment_overlaps(beta,frame,momenta,shape):
    """Position jets use actual placed beta and physical p, no p derivatives."""
    result = {}
    for axes in MOMENTS:
        monomial = np.prod(momenta[:,axes],axis=1) if axes else np.ones(len(momenta))
        result[axes] = zgemm(1.,beta,monomial[:,None]*frame,trans_a=2).reshape(shape)
    return result


def product_third(left,right,active_d):
    """Trace of all eight third product legs for independent complex dyads.

    left/right are B†p^axes L/R overlaps [atom,channel,column], where
    rho=L R†. The D below is UPF's expanded NC coupling, not state D.
    """
    out = np.zeros((len(left[()]),3,3,3),complex)
    dl = {axes: np.einsum('ij,ajb->aib',active_d,value,optimize=True)
          for axes,value in left.items()}
    for axes in TRIPLES:
        value = np.zeros(len(out),complex)
        for mask in range(8):
            bra_axes = tuple(axes[leg] for leg in range(3) if mask & (1<<leg))
            ket_axes = tuple(axes[leg] for leg in range(3) if not mask & (1<<leg))
            phase = (-1j)**len(bra_axes)*(1j)**len(ket_axes)
            value += phase*np.einsum('aib,aib->a',right[bra_axes].conj(),dl[ket_axes],optimize=True)
        for a,b,c in set(permutations(axes)):
            out[:,a,b,c] = value
    return out


def local_third(density,coupling):
    """Omega sum_G n(-G) (+i G_a G_b G_c) vloc(G) exp(-iG.tau)."""
    ionic,fft = coupling.ionic,coupling.factors.fft
    if np.shape(density) != (fft.nnr,):
        raise ValueError('original full physical total density required')
    transformed = np.fft.ifftn(np.asarray(density).reshape(fft.physical,order='F'))
    miller = ionic.inputs.density.miller
    sampled = transformed[tuple(miller[i]%fft.physical[i] for i in range(3))]
    g = ionic.g; ref = ionic.inputs.reference
    atoms = tuple(getattr(coupling,'cubic_atom_indices',range(ref.nat)))
    out = np.zeros((len(atoms),3,3,3),complex)
    moduli = np.linalg.norm(g,axis=1)
    for index,atom in enumerate(atoms):
        name = ref.atom_names[atom]
        if name not in ionic.vloc:
            ionic.vloc[name] = ionic.inputs.tables[name].vloc(moduli)
        weight = sampled*ionic.vloc[name]*np.exp(-1j*g@ref.positions_bohr[:,atom])
        out[index] = 1j*fft.omega*np.einsum('g,ga,gb,gc->abc',weight,g,g,g,optimize=True)
    return out


def full_column_projector_third(psi,coefficient,coupling):
    """Associate B†p^axes(Psi D) as (B†p^axes Psi)D, without PW replay.

    Every raw column and every complex D entry remains present. Full-column
    moment GEMMs are reused on both product legs; workspace is20*P*B rather
    than a full Psi D or a PW-square operator.
    """
    coverage = tuple(getattr(coupling,'cubic_atom_indices',range(coupling.ionic.inputs.reference.nat)))
    rowmap = {atom:index for index,atom in enumerate(coverage)}
    out = np.zeros((len(coverage),3,3,3),complex)
    for species in coupling.factors.species:
        atoms = np.asarray(species.atom_indices,int)-1
        shape = (len(atoms),len(species.active_d),psi.shape[1])
        right = (compressed_overlaps(species,psi,coupling.ionic.p,
                    coupling.cubic_radials[species.name],MOMENTS)
                 if hasattr(coupling,'cubic_radials')
                 else moment_overlaps(species.beta,psi,coupling.ionic.p,shape))
        left = {axes:zgemm(1.,value.reshape(-1,psi.shape[1]),coefficient).reshape(shape)
                for axes,value in right.items()}
        out[[rowmap[int(atom)] for atom in atoms]] = product_third(left,right,species.active_d)
    return out


def full_state_third(psi,coefficient,coupling):
    """One full original sector, both local and nonlocal, column-blocked.

    Returns unweighted local/nonlocal[declared_atoms,3,3,3]. Atom coverage is
    explicit, never zero-filled; all original raw columns and complex D survive.
    """
    if (psi.ndim != 2 or psi.shape[0] != coupling.npw
            or coefficient.shape != (psi.shape[1],psi.shape[1])):
        raise ValueError('complete raw complex frame and square state coefficient required')
    fft = coupling.factors.fft
    density = np.zeros(fft.nnr,complex)
    nonlocal_term = full_column_projector_third(psi,coefficient,coupling)
    for start in range(0,psi.shape[1],coupling.block_size):
        stop = min(start+coupling.block_size,psi.shape[1])
        left = zgemm(1.,psi,coefficient[:,start:stop])
        right = psi[:,start:stop]
        density += np.sum(position_fft(fft,left)*position_fft(fft,right).conj(),axis=1)/fft.omega
    return local_third(density,coupling),nonlocal_term
