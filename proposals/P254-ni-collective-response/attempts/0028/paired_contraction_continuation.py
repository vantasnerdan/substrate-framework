"""Retain an actual passed forward witness across only the dual factor fix."""
import ast
from copy import deepcopy
from pathlib import Path

import sympy as sp

from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_material_state import authenticate, file_binding

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def exact_dual_factorization(check):
    nt, ns, rank = sp.symbols('nt ns rank', integer=True, positive=True)
    target = sp.MatrixSymbol('target', nt, rank)
    response = sp.MatrixSymbol('response', ns, rank)
    coefficient = sp.MatrixSymbol('coefficient', rank, rank)
    check('SymPy_arbitrary_rank_complex_dual_factor_identity',
          sp.simplify((-sp.I*target*coefficient)*sp.adjoint(response)
                      - target*sp.adjoint(sp.I*response*sp.adjoint(coefficient))) == sp.ZeroMatrix(nt, ns))
    p = sp.Matrix([[1+sp.I, 2], [-sp.I, 3-sp.I], [2, 1+2*sp.I]])
    r = sp.Matrix([[2, sp.I], [1-sp.I, 3], [-1, 2+sp.I], [sp.I, 1]])
    d = sp.Matrix([[2+sp.I, 1-2*sp.I], [3+sp.I, -1]])
    old_dyad = (-sp.I*p*d)*r.H
    new_dyad = p*(sp.I*r*d.H).H
    check('SymPy_dual_dyad_factor_transfer_without_Hermitian_assumption',
          all(sp.expand(z) == 0 for z in old_dyad-new_dyad))
    check('SymPy_missing_factor_i_mutation_rejected',
          any(sp.expand(z) != 0 for z in old_dyad-p*(r*d.H).H))
    check('SymPy_missing_D_adjoint_mutation_rejected',
          any(sp.expand(z) != 0 for z in old_dyad-p*(sp.I*r*d).H))
    a = sp.Matrix([[1, sp.I, 2, -1], [2-sp.I, 1, sp.I, 3], [-sp.I, 2, 1, 2+sp.I]])
    check('SymPy_weighted_moment_leg_is_not_original_adjoint_drive_frame',
          any(sp.expand(z) != 0 for z in a.H*p-a.H*(-sp.I*p*d)))


def _dump(tree):
    return ast.dump(tree, include_attributes=False)


def _only_dual_factor_changed(before, after):
    old_statement = ast.parse('left2, right2 = -1j*pt @ t.d, x.bra[index].conj().T').body[0]
    new_statement = ast.parse('left2, right2 = pt, 1j*x.bra[index].conj().T @ t.d.conj().T').body[0]
    old, new = ast.parse(before), ast.parse(after)
    class Restore(ast.NodeTransformer):
        count = 0
        def visit_Assign(self, node):
            if _dump(node) == _dump(new_statement):
                self.count += 1
                return deepcopy(old_statement)
            return self.generic_visit(node)
    restore = Restore()
    new = restore.visit(new)
    return restore.count == 1 and _dump(old) == _dump(new)


def _only_proof_inventory_changed(before, after):
    old, new = ast.parse(before), ast.parse(after)
    expected = _dump(ast.parse("HERE/'paired_contraction_continuation.py'", mode='eval').body)
    class Remove(ast.NodeTransformer):
        count = 0
        def visit_List(self, node):
            node = self.generic_visit(node)
            retained = [item for item in node.elts if _dump(item) != expected]
            self.count += len(node.elts)-len(retained)
            node.elts = retained
            return node
    remove = Remove()
    new = remove.visit(new)
    return remove.count == 1 and _dump(old) == _dump(new)


def retain_forward(args, sources, initialization, history, check):
    binding = file_binding(args.retained_forward)
    registration = file_binding(args.previous_registration)
    if binding['sha256'] != args.retained_forward_sha256 or registration['sha256'] != args.previous_registration_sha256:
        raise ValueError('explicit retained forward/registration digest differs')
    previous, registered = bound_json(binding), bound_json(registration)
    prior = {source['path']: source for source in previous['source_bindings']}
    current = {source['path']: source for source in sources}
    edge_path = str(ROOT/'src/substrate_framework/partitioned_q_response.py')
    driver_path = str(HERE/'run_joint_q_material.py')
    witness_path = str(HERE/'run_paired_contraction_consumer_check.py')
    helper_path = str(Path(__file__).resolve())
    check('retained_forward_source_delta_is_only_dual_and_proof_inventory',
          set(current)-set(prior) == {helper_path} and set(prior).issubset(current)
          and {p for p in prior if prior[p]['sha256'] != current[p]['sha256']}
          == {edge_path, driver_path, witness_path})
    for path, name, comparator in (
        (edge_path, 'partitioned_q_response.pre-dual-factor-fix-001.py', _only_dual_factor_changed),
        (driver_path, 'run_joint_q_material.pre-dual-factor-fix-001.py', _only_proof_inventory_changed)):
        snapshot = file_binding(HERE/'source-history'/name)
        check(name+'_exact_prior_executable_bytes', snapshot['sha256'] == prior[path]['sha256'])
        check(name+'_only_declared_nonforward_change', comparator(authenticate(snapshot).read_bytes(), Path(path).read_bytes()))
    witness = file_binding(HERE/'source-history/run_paired_contraction_consumer_check.pre-dual-factor-fix-001.py')
    check('original_forward_probe_source_preserved_exactly', witness['sha256'] == prior[witness_path]['sha256'])
    def same_binding(a, b):
        return a['path'] == b['path'] and a['sha256'] == b['sha256'] and str(a['bytes']) == str(b['bytes'])
    command = registered['command']
    check('retained_forward_has_identical_original_material_inputs_and_target',
          same_binding(initialization, registered['initialization']) and same_binding(history, registered['history'])
          and command[command.index('--case')+1] == args.case
          and int(command[command.index('--edge-index')+1]) == args.edge_index
          and float(command[command.index('--absolute-target')+1]) == args.absolute_target)
    check('retained_actual_forward_is_complete_not_a_dual_or_full_admission',
          previous['status'] == 'partial-until-both-directions'
          and all(item['passed'] is True for item in previous['checks'])
          and all(any(item['name'] == 'forward_actual_full_edge_'+name and item['passed'] is True
                      for item in previous['checks']) for name in ('ket', 'bra', 'n', 'g', 'ell'))
          and [item['direction'] for item in previous['timings']] == ['forward'])
    return dict(progress=binding, original_registration=registration, original_probe_source=witness,
        timings=previous['timings'], checks=previous['checks'], rerun=False,
        scope='Actual original forward witness retained because the entire executable outside the one dual factor assignment and proof inventory is unchanged. The corrected dual must execute independently. Old failed result is not relabelled.')
