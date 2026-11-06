"""Source CODATA2006 conversion against an independent symbolic SI route."""
import ast
from pathlib import Path
import numpy as np
import sympy as sp
from substrate_framework.material_primary_boundary import atomic_units_2006


def run(ledger, context):
    source = Path(context['supplement']['ase_root'])/'units.py'
    tree = ast.parse(source.read_text())
    table = next(n.value for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='CODATA' for t in n.targets))
    section = next(v for k,v in zip(table.keys,table.values) if ast.literal_eval(k)=='2006')
    constants = {ast.literal_eval(k):sp.sympify(ast.unparse(v),locals={'pi':sp.pi},rational=True) for k,v in zip(section.keys,section.values)}
    a,ry = atomic_units_2006(*(constants[k] for k in ('_c','_mu0','_hplanck','_e','_me')))
    # Independent route: fine structure alpha and rest energy.
    c,h,e,me,mu = (constants[k] for k in ('_c','_hplanck','_e','_me','_mu0'))
    alpha = mu*c*e**2/(2*h)
    other_a = h/(2*sp.pi*me*c*alpha)*10**10
    other_ry = me*c**2*alpha**2/(2*e)
    ledger.check('SI Coulomb and fine-structure conversions coincide symbolically', sp.simplify(a-other_a)==0 and sp.simplify(ry-other_ry)==0)
    from ase.units import create_units
    ase = create_units('2006')
    af,rf = float(a),float(ry)
    allowance = 32*np.finfo(float).eps
    ledger.mutation_sensitive('source-specific conversion differs from Hartree and wrong Bohr convention',
        lambda x: abs(x[0]/ase['Bohr']-1)<allowance and abs(x[1]/ase['Ry']-1)<allowance,
        (af,rf), [(af,2*rf),(af*10,rf),(af,create_units('2014')['Ry'])])
    ledger.check('ASE espresso source explicitly selects CODATA2006', "units = create_units('2006')" in (source.parent/'io/espresso.py').read_text())
    context['units']={'bohr_angstrom':af,'rydberg_ev':rf,'force_factor':rf/af,'elementary_charge_si':float(e)}
    return {**context['units'],'constants':{k:str(v) for k,v in constants.items()},
            'convention':'ASE espresso CODATA2006, DOI 10.1103/RevModPhys.80.633',
            'symbolic_bohr_angstrom':str(a),'symbolic_rydberg_ev':str(ry),
            'roundoff_relative_allowance':allowance,
            'scope':'Conversion of published stdout decimals, not native numerical convergence error.'}
