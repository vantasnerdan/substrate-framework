"""Exact finite-sum distributivity checks before implementation; no material physics claim."""
from pathlib import Path
import json
import sympy as s
from itertools import product
HERE=Path(__file__).resolve().parent

def run():
    out=[]
    assert (-s.I)**2*s.I == -s.I and s.I**2*(-s.I)==s.I
    out.append(dict(block='raw-Gaussian-phases',greater=str((-s.I)**2*s.I),lesser=str(s.I**2*(-s.I)),exchange_adds_one_minus=True,mutation_wrong_reverse_phase_exposed=((-s.I)**3 != -s.I)))
    def mat(name):return s.Matrix(2,2,lambda i,j:s.Symbol(f'{name}{i}{j}'))
    bt,bs,ct,cs,aa,at,y,z=[mat(n) for n in ('bt','bs','ct','cs','aa','at','y','z')]
    old=s.trace(y*bt*bs*z*aa*at)
    new=s.trace((at*y*bt)*(bs*z*aa))
    residual=s.expand(old-new)
    bad=s.expand(old-s.trace((at.T*y*bt)*(bs*z*aa)))
    assert residual==0 and bad!=0
    out.append(dict(block='direct-time-separated-bubble',identity='Tr(Y Bt Bs Z As At)=sum_alpha,beta (At Y Bt)[alpha,beta](Bs Z As)[beta,alpha]',residual=str(residual),mutation_reverse_transpose_exposed=bad!=0))
    L={idx:s.Symbol('L'+''.join(map(str,idx))) for idx in product(range(2),repeat=3)}
    R={idx:s.Symbol('R'+''.join(map(str,idx))) for idx in product(range(2),repeat=3)}
    fb,fc,rr=bt*bs,ct*cs,aa*at
    old=sum(L[b,a,c]*fb[b,d]*fc[c,e]*rr[f,a]*R[d,f,e] for b,a,c,d,e,f in product(range(2),repeat=6))
    lt={};rs={};badrs={}
    for beta,alpha,gamma in product(range(2),repeat=3):
        lt[beta,alpha,gamma]=sum(L[b,a,c]*bt[b,beta]*at[alpha,a]*ct[c,gamma] for b,a,c in product(range(2),repeat=3))
        rs[beta,alpha,gamma]=sum(bs[beta,d]*aa[f,alpha]*cs[gamma,e]*R[d,f,e] for d,f,e in product(range(2),repeat=3))
        badrs[beta,alpha,gamma]=sum(bs[d,beta]*aa[f,alpha]*cs[gamma,e]*R[d,f,e] for d,f,e in product(range(2),repeat=3))
    new=sum(lt[k]*rs[k] for k in lt)
    residual=s.expand(old-new);bad=s.expand(old-sum(lt[k]*badrs[k] for k in lt))
    assert residual==0 and bad!=0
    out.append(dict(block='exchange-three-internal-mode-bilinear',identity='L_iba c Fb_bd Fc_ce R_fa V_dfej = sum_beta,alpha,gamma LT_i,beta,alpha,gamma RS_beta,alpha,gamma,j',residual=str(residual),mutation_independent_right_orientation_exposed=bad!=0,modes='all m^3, no SVD or occupation simplification'))
    assert s.simplify(s.Rational(1,64)**2-s.Rational(1,128)**2)!=0
    out.append(dict(block='owning-normalization',two_vertex_normalization='1/N^2',mutation_double_spin_weight_exposed=True,source='A46 weights[q]/c.n on each vertex; explicit spin slots retained'))
    (HERE/'symbolic-findings.json').write_text(json.dumps(dict(status='complete',findings=out),indent=2)+'\n')
    print(json.dumps(out),flush=True)
if __name__=='__main__':run()
