"""Exact oriented-coupling/tau fixture for the discrete QE force-lattice convention."""

from copy import deepcopy
import numpy as np
import sympy as sp

from substrate_framework import material_force_lattice as lattice


def zero(matrix):
    return all(sp.simplify(x) == 0 for x in matrix)


def force_lattice_block(ledger):
    # Real, reciprocal finite Hessian with nonsymmetric oriented Cartesian bonds.
    a = sp.Matrix([[2,1,0],[0,3,2],[1,0,4]])
    k = a + a.T
    c0 = sp.diag(k,k)
    c1 = sp.zeros(6)
    c1[:3,3:] = -a
    c1[3:,:3] = -a
    c2 = c1.T
    realspace = [c0,c1,c2]
    omega = -sp.Rational(1,2) + sp.I*sp.sqrt(3)/2
    reciprocal = [sum((realspace[r]*omega**(-q*r) for r in range(3)),sp.zeros(6)).applyfunc(sp.simplify)
                  for q in range(3)]
    recovered = [sum((reciprocal[q]*omega**(q*r)/3 for q in range(3)),sp.zeros(6)).applyfunc(sp.simplify)
                 for r in range(3)]
    ledger.check("exact nontrivial discrete Fourier reconstruction", all(zero(x-y) for x,y in zip(recovered,realspace)))
    ledger.check("exact reciprocal oriented bonds", all(zero(realspace[r]-realspace[-r % 3].T) for r in range(3)))
    translations = sp.eye(3).col_join(sp.eye(3))
    ledger.check("exact force-lattice translations", zero(reciprocal[0]*translations))
    ledger.check("exact Hermitian q matrices", all(zero(q-q.conjugate().T) for q in reciprocal))
    wrong_sign = [sum((reciprocal[q]*omega**(-q*r)/3 for q in range(3)),sp.zeros(6)).applyfunc(sp.simplify)
                  for r in range(3)]
    ledger.check("Fourier sign mutation exposes oriented bond transpose", not zero(wrong_sign[1]-realspace[1]))
    tau = [sp.zeros(3,1),sp.Matrix([sp.Rational(1,5),sp.Rational(2,7),sp.Rational(1,3)])]
    separation = sp.Matrix([1,0,0])+tau[0]-tau[1]
    ledger.check("oriented geometry includes tau difference", separation != sp.Matrix([1,0,0])+tau[1]-tau[0])
    # At q=(1/3,0,0), incorrectly attaching tau to the native phase multiplies
    # an off-diagonal atom block by exp(±2pi i/15), not unity.
    theta = 2*sp.pi/15
    ledger.check("incorrect atom-phase mutation is nonidentity", sp.simplify(1-sp.cos(theta)).is_positive is True)
    expected_r = np.stack([np.asarray(x,dtype=complex).reshape(2,3,2,3) for x in realspace])[:,None,None]
    expected_q = np.stack([np.asarray(x,dtype=complex).reshape(2,3,2,3) for x in reciprocal])[:,None,None]
    result_r = lattice.q_grid_to_force_lattice(expected_q)
    tolerance = 128*np.finfo(float).eps*np.linalg.norm(expected_r)
    ledger.check("numerical transform matches independent exact sums", np.linalg.norm(result_r-expected_r) <= tolerance)
    ledger.check("numerical inverse reconstructs exact q grid", np.linalg.norm(lattice.force_lattice_to_q_grid(result_r)-expected_q) <= tolerance)
    tau_array = np.array([[0,0,0],[1/5,2/7,1/3]])
    wrong_q = expected_q.copy()
    for q in range(3):
        phase = np.exp(-2j*np.pi*q/3*(tau_array[:,0,None]-tau_array[None,:,0]))
        wrong_q[q,0,0] *= phase[:,None,:,None]
    ledger.check("tau-phase mutation fails actual force-lattice reconstruction", np.linalg.norm(lattice.q_grid_to_force_lattice(wrong_q)-expected_r) > tolerance)
    residuals = lattice.lattice_residuals(result_r)
    ledger.check("force lattice retains real reciprocity and translation", max(residuals.values()) <= tolerance)
    corrupt = result_r.copy()
    corrupt[1,0,0,0,0,1,1] += 1j
    ledger.check("complex reciprocity mutation exposed", lattice.lattice_residuals(corrupt)["reciprocity_frobenius"] > tolerance)
    geometry = {"number_of_atoms":2,"at_cartesian_by_lattice_vector_alat":np.eye(3).tolist()}
    members = [{"q_cartesian_2pi_alat":[q/3,0,0],
                "force_matrix_ry_per_bohr2":{"real":expected_q[q,0,0].real.tolist(),"imag":expected_q[q,0,0].imag.tolist()}}
               for q in range(3)]
    record = {"geometry":geometry,"q_star":members}
    completed = lattice.complete_native_q_grid([record],[3,1,1],grid_atol=1e-12)
    ledger.check("complete q-star expansion preserves every complex entry", np.array_equal(completed,expected_q))
    for name,change in (("missing",members[:-1]),("duplicate",members+[members[0]])):
        bad = deepcopy(record)
        bad["q_star"] = change
        try:
            lattice.complete_native_q_grid([bad],[3,1,1],grid_atol=1e-12)
        except ValueError:
            rejected = True
        else:
            rejected = False
        ledger.check(f"{name} q-star mutation rejected", rejected)
    return {"tau_fractional": [[str(v) for v in t] for t in tau],
            "oriented_cartesian_coupling":str(a),"numeric_residuals":residuals,
            "scope":"Exact synthetic 3x1x1 periodic lattice only; no infinite-range pair law or off-grid interpolation."}
