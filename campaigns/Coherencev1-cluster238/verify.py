"""Verify conditional contact-shift references without a material conclusion."""
import argparse
import json
from pathlib import Path
import time
import sympy as sp


def verify():
    C, R, Rref, n = sp.symbols('C R Rref n',positive=True)
    initial, united = C/R, 2*C/R
    change = (initial-united)-(C/Rref-2*C/Rref)
    expected = C*(1/Rref-1/R)
    assert sp.simplify(change-expected) == 0
    assert sp.simplify(sp.diff(expected,R)-C/R**2) == 0
    U = n**sp.Rational(1,6)
    assert sp.simplify(n*sp.diff(U,n)/U-sp.Rational(1,6)) == 0
    assert sp.simplify(expected.subs(R,Rref)) == 0
    # Wrong positive C/R interpretation of a nearby positive spectator is rejected.
    assert sp.simplify((C/R-C/Rref)-expected) != 0
    return dict(status='passed',checks=['balanced initial/united spectator energy',
         'closer positive spectator increases barrier','matched reference cancellation',
         'theoretical TF density exponent1/6'],
         positive_spectator_screening_mutation_rejected=True,
         scope='two separate conditional closures; not a material dielectric calculation')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    started=time.process_time()
    result=verify()
    result['process_cpu_s']=time.process_time()-started
    Path(args.output).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result))
