"""Exact NC position moments, sharing polynomial/radial work without interpolation."""
from functools import lru_cache
import numpy as np
from scipy.fft import ifftn
from scipy.linalg.blas import zgemm
from substrate_framework.qe_bloch_frame import BlochBasis
from substrate_framework.qe_real_bloch_grid import RealBlochBasis
from substrate_framework.qe_material_grid import PhysicalRhoFFT
from substrate_framework.qe_material_hamiltonian import PWFactors,SpeciesProjectors
from substrate_framework.qe_material_field_derivatives import UPFRadial,solid_harmonics
from substrate_framework.qe_nc_projectors import place_projectors
from substrate_framework.material_physical_operators import REPRESENTATION


@lru_cache(maxsize=1)
def solid_polynomials():
    """Extract the exact degree<=2 polynomial from the established angular API."""
    zero=np.zeros((1,3))
    evaluate=lambda p:solid_harmonics(np.asarray(p).reshape(1,3),np.zeros(3),np.zeros(3))[0][0]
    constant=evaluate(zero)
    coefficients={(0,0,0):constant}
    eye=np.eye(3)
    for a in range(3):
        exponent=tuple(int(i==a) for i in range(3))
        coefficients[exponent]=(evaluate(eye[a])-evaluate(-eye[a]))/2
        coefficients[tuple(2*x for x in exponent)]=(evaluate(eye[a])+evaluate(-eye[a]))/2-constant
        for b in range(a):
            coefficients[tuple(int(i in (a,b)) for i in range(3))]=evaluate(eye[a]+eye[b])-evaluate(eye[a])-evaluate(eye[b])+constant
    return tuple({power:float(value[lm]) for power,value in coefficients.items() if value[lm]!=0.} for lm in range(9))


def position_factors(inputs,basis,electronic_atoms=None):
    """The established placed physical beta formula, deduplicating exact s only.

    No interpolation, rounding of s, H/F/J or radial derivative is introduced.
    The original real-basis/UPF identity checks remain compulsory.
    """
    if not isinstance(basis,(BlochBasis,RealBlochBasis)):
        raise ValueError('checked original physical basis required')
    ref=inputs.reference; basis.validate(ref)
    selected=set(range(ref.nat) if electronic_atoms is None else electronic_atoms)
    real=isinstance(basis,RealBlochBasis)
    if real:
        actual={(name,table.upf.sha256) for name,table in inputs.tables.items()}
        expected={(s.name,next(x['sha256'] for x in basis.scope.upf_sources if x['path'].endswith('/'+s.pseudo_file))) for s in ref.species}
        if actual!=expected:
            raise ValueError('physical factor UPFs differ from admitted real scope')
    volume=abs(float(np.linalg.det(ref.cell_bohr)))
    fft=PhysicalRhoFFT.create(ref.fft_smooth,basis.miller,basis.reciprocal_bohr_inv,volume)
    p=(basis.reciprocal_bohr_inv@basis.miller+basis.k_cart_bohr_inv[:,None]).T
    squared=np.einsum('pi,pi->p',p,p)
    unique,inverse=np.unique(squared,return_inverse=True)
    solid=solid_harmonics(p,np.zeros(3),np.zeros(3))[0]
    groups=[]; compact={}
    for species in ref.species:
        table=inputs.tables[species.name]; upf=table.upf
        if table.omega_bohr3!=volume:
            raise ValueError('original radial volume differs')
        radial=UPFRadial(table).evaluate(unique)[0][inverse]
        center=solid[:,upf.channel_lm]*radial[:,upf.channel_radial]*((-1j)**upf.channel_l)[None,:]
        atoms=tuple(i for i,name in enumerate(ref.atom_names) if name==species.name and i in selected)
        if not atoms:
            continue
        beta=np.empty((len(p),len(atoms)*upf.nh),complex,order='F')
        for slot,atom in enumerate(atoms):
            beta[:,slot*upf.nh:(slot+1)*upf.nh]=place_projectors(center,p,ref.positions_bohr[:,atom])
        groups.append(SpeciesProjectors(species.name,tuple(a+1 for a in atoms),beta,upf.expanded_d))
        compact[species.name]=(radial,upf,ref.positions_bohr[:,atoms])
    provenance=dict(kind='real-collinear-orbital-conjugation-basis' if real else 'unitary-reconstructed-basis',identity=basis.identity,
        source=dict(basis.slot.source_binding),xml_source=dict(basis.reference_source),ispin=basis.slot.ispin,mesh_index=basis.slot.mesh_index,
        source_ik=basis.slot.source_ik,operation=basis.slot.operation.index,gauge_policy=basis.slot.gauge_policy,representation=REPRESENTATION)
    if real:
        provenance.update(scope=basis.scope.identity,orbital_conjugation=basis.slot.orbital_conjugation)
    return PWFactors(fft,squared,tuple(groups),basis.k_cart_bohr_inv,None,inputs.sources+(provenance,)),compact


def compressed_overlaps(species,frame,p,compact,moments):
    radial,upf,positions=compact
    nat=len(species.atom_indices); bands=frame.shape[1]
    out={axes:np.zeros((nat,len(species.active_d),bands),complex) for axes in moments}
    powers={axes:tuple(axes.count(i) for i in range(3)) for axes in moments}
    polynomials=solid_polynomials()
    atomic_phase=np.exp(-1j*p@positions)
    for radial_index in range(upf.nbeta):
        channels=np.flatnonzero(upf.channel_radial==radial_index)
        required=sorted({tuple(a+b for a,b in zip(power,powers[axes])) for channel in channels
            for power in polynomials[int(upf.channel_lm[channel])] for axes in moments})
        lookup={power:index for index,power in enumerate(required)}
        bank=np.empty((len(p),len(required)*nat),complex,order='F')
        base=radial[:,radial_index,None]*atomic_phase
        for index,power in enumerate(required):
            bank[:,index*nat:(index+1)*nat]=base*np.prod(p**np.asarray(power),axis=1)[:,None]
        projected=zgemm(1.,bank,frame,trans_a=2).reshape(len(required),nat,bands)
        for channel in channels:
            phase=(1j)**int(upf.channel_l[channel])
            for axes in moments:
                for power,coefficient in polynomials[int(upf.channel_lm[channel])].items():
                    exponent=tuple(a+b for a,b in zip(power,powers[axes]))
                    out[axes][:,channel]+=phase*coefficient*projected[lookup[exponent]]
    return out


def position_fft(fft,coefficients):
    """Same explicit original scatter/normalization, single-thread in-place FFT."""
    spectrum=np.zeros((fft.nnr,*coefficients.shape[1:]),complex)
    spectrum[fft.nl]=coefficients
    field=ifftn(fft.physical_view(spectrum),axes=(0,1,2),norm='forward',workers=1,overwrite_x=True)
    return field.reshape(spectrum.shape,order='F')
