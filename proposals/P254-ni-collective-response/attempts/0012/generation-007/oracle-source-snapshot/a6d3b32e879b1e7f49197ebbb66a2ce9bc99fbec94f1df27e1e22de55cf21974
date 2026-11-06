"""Check source CODATA2006 by independent fine-structure and executable ASE routes."""
import ast
import numpy as np
import sympy as sp
from substrate_framework.material_primary_boundary import atomic_units_2006
import source_pi


def source_constants(text,year):
    tree=ast.parse(text)
    table=next(n.value for n in tree.body if isinstance(n,ast.Assign) and
               any(isinstance(t,ast.Name) and t.id=='CODATA' for t in n.targets))
    section=next(v for k,v in zip(table.keys,table.values) if ast.literal_eval(k)==year)
    return {ast.literal_eval(k):sp.sympify(ast.get_source_segment(text,v),
            locals={'pi':sp.pi},rational=True) for k,v in zip(section.keys,section.values)}


def espresso_functions(text,path,create_units):
    tree=ast.parse(text)
    selected=[n for n in tree.body if (isinstance(n,ast.Assign) and
              any(isinstance(t,ast.Name) and t.id=='units' for t in n.targets)) or
              (isinstance(n,ast.FunctionDef) and n.name in ('_get_forces','_get_energy'))]
    namespace={'create_units':create_units,'np':np}
    exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),namespace)
    return namespace


def source_functions(sources,units_ref,espresso_ref,name):
    text=sources.text(units_ref['path'],units_ref['sha256'])
    module=sources.module(units_ref['path'],units_ref['sha256'],name)
    espresso_text=sources.text(espresso_ref['path'],espresso_ref['sha256'])
    espresso=espresso_functions(espresso_text,espresso_ref['path'],module.create_units)
    return text,espresso


def run(ledger,context):
    registered=context['inventory']['raw_registered_sources']
    units_ref=next(r for r in registered if r['path'].endswith('/ase/units.py'))
    espresso_ref=next(r for r in registered if r['path'].endswith('/ase/io/espresso.py'))
    text,espresso=source_functions(context['sources'],units_ref,espresso_ref,'p254_registered_ase_units')
    constants=source_constants(text,'2006')
    a,ry=atomic_units_2006(*(constants[k] for k in ('_c','_mu0','_hplanck','_e','_me')))
    c,h,e,me,mu=(constants[k] for k in ('_c','_hplanck','_e','_me','_mu0'))
    alpha=mu*c*e**2/(2*h)
    independent_a=h/(2*sp.pi*me*c*alpha)*10**10
    independent_ry=me*c**2*alpha**2/(2*e)
    ledger.check('source SI Coulomb and independent fine-structure routes agree exactly',
                 sp.simplify(a-independent_a)==0 and sp.simplify(ry-independent_ry)==0)
    certified=source_pi.run(ledger,a,ry)
    context['force_factor_bounds']=certified['force_bounds']
    from check_source_consumption import run as check_consumption
    source_proof=check_consumption(ledger,context,units_ref,espresso_ref)
    numeric=(espresso['units']['Bohr'],espresso['units']['Ry'])
    allowance=32*np.finfo(float).eps
    predicate=lambda pair: all(abs(float(exact)/observed-1)<allowance
                               for exact,observed in zip(pair,numeric))
    altered=source_constants(text,'2014')
    other=atomic_units_2006(*(altered[k] for k in ('_c','_mu0','_hplanck','_e','_me')))
    ledger.mutation_sensitive('numerical actual source-defined espresso convention',predicate,
                              (a,ry),[(a,2*ry),(a*10,ry),other])
    context['units']={'bohr':a,'ry':ry,'force':sp.cancel(ry/a)}
    context['espresso']=espresso
    return {'constants_exact':{k:str(v) for k,v in constants.items()},
            'verified_source_consumer_proof':source_proof,
            'actual_source_bindings':{role:context['sources'].record(ref['path'],ref['sha256'])
                                     for role,ref in [('units',units_ref),('espresso',espresso_ref)]},
            'certified_force_factor_bounds_exact':[str(v) for v in certified['force_bounds']],
            'certified_pi_bounds_exact':[str(v) for v in certified['pi_bounds']],
            'certified_conversion_proof':certified['proof'],
            'bohr_angstrom_exact':str(a),'rydberg_ev_exact':str(ry),'force_factor_exact':str(sp.cancel(ry/a)),
            'bohr_angstrom':float(a),'rydberg_ev':float(ry),'force_factor':float(ry/a),
            'actual_source_espresso_units':{k:numeric[i] for i,k in enumerate(('Bohr','Ry'))},
            'route':'fine structure alpha/rest energy, independent of reused Coulomb route; source AST executed, no wording assertion',
            'floating_relative_crosscheck_allowance':allowance,
            'scope':'Source-defined CODATA2006 convention; constants are fixed inputs, not physical uncertainty bounds.'}
