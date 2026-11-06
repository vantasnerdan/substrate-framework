"""Compare independently decoded native XML spectra with every raw q-star matrix."""

from pathlib import Path
import json
import xml.etree.ElementTree as ET
import numpy as np

from substrate_framework import material_modes as modes
from substrate_framework import material_force_lattice as lattice


def xml_numbers(node):
    return np.array([float(x.replace("D","E")) for x in node.text.split()])


def xml_force(node, nat):
    """Independent entrywise decode: Fortran Cartesian column index is outermost."""
    force = np.empty((nat,3,nat,3),dtype=complex)
    for i in range(nat):
        for j in range(nat):
            pairs = xml_numbers(node.find(f"PHI.{i+1}.{j+1}")).reshape(9,2)
            for b in range(3):
                for a in range(3):
                    force[i,a,j,b] = complex(*pairs[3*b+a])
    return force


def spectral_budget(matrix, reference_lambda, operation_scale):
    """Engineering float64 allowance, including the full-grid transform scale.

    The sum of sampled operator norms bounds the scale entering Fourier sums;
    the small residual Gamma matrix alone would miss cancellation roundoff.
    This is not a certified LAPACK enclosure or a material convergence error.
    """
    n = len(reference_lambda)
    scale = max(np.linalg.norm(matrix,2),np.max(np.abs(reference_lambda)),operation_scale)
    return float(256*n*np.finfo(float).eps*scale)


def frequency_budget(values, reference, lambda_budget):
    """Propagate signed-sqrt absolute errors, avoiding relative acoustic tests."""
    same_sign = values*reference > 0
    denominator = np.sqrt(np.abs(values)) + np.sqrt(np.abs(reference))
    budget = np.full(values.shape,np.sqrt(2*lambda_budget))
    budget[same_sign] = lambda_budget/denominator[same_sign]
    return budget*modes.RY_TO_THZ


def frequency_interval_qe(frequency_thz, half_width_thz):
    """Invert the signed square-root monotonically, including intervals across zero."""
    if not np.isfinite(half_width_thz) or half_width_thz < 0:
        raise ValueError("frequency rounding width must be finite and nonnegative")
    low = np.asarray(frequency_thz) - half_width_thz
    high = np.asarray(frequency_thz) + half_width_thz
    return (np.sign(low)*(low/modes.RY_TO_THZ)**2,
            np.sign(high)*(high/modes.RY_TO_THZ)**2)


def degeneracy_clusters(values, absolute_budget):
    clusters = [[0]]
    for index in range(1,len(values)):
        if values[index]-values[index-1] <= 2*absolute_budget:
            clusters[-1].append(index)
        else:
            clusters.append([index])
    return clusters


def native_pilot_block(ledger, result_path, native_directory, replay):
    native_directory = Path(native_directory)
    envelope = json.loads(Path(result_path).read_text())
    ledger.check("native pilot published complete", envelope["status"] == "complete")
    data = envelope["data"]
    grid_shape = data["q_manifest"]["grid"]
    grid = lattice.complete_native_q_grid(data["dynamical_matrices"],grid_shape,grid_atol=1e-12)
    masses = modes.atom_masses_amu(data["dynamical_matrices"][0]["geometry"])
    operation_scale = sum(np.linalg.norm(modes.mass_weighted_operator(grid[index],masses,units="qe_ry"),2)
                          for index in np.ndindex(*grid_shape))
    replay_rows = {tuple(row["grid_index"]):row for row in replay["rows"]}
    ledger.check("native replay covers each sampled q exactly once",
                 set(replay_rows) == set(np.ndindex(*grid_shape)) and len(replay_rows) == len(replay["rows"]))
    cm_1_to_thz = modes.C_SI/1e10
    # matdyn:779 is (6f10.4), not an empirical fitted tolerance.
    half_width_thz = (0.5*10**-4)*cm_1_to_thz
    rows = []
    for record in data["dynamical_matrices"]:
        native = ET.parse(native_directory/Path(record["artifact"]).name).getroot()
        nat = record["geometry"]["number_of_atoms"]
        at = np.asarray(record["geometry"]["at_cartesian_by_lattice_vector_alat"])
        ph_reported = np.sort([xml_numbers(native.find(f"FREQUENCIES_THZ_CMM1/OMEGA.{i+1}"))[0]
                               for i in range(3*nat)])
        xml_stars = [node for node in native if node.tag.startswith("DYNAMICAL_MAT_.")]
        ledger.check(f"q{record['q_index']} complete XML star count", len(xml_stars) == len(record["q_star"]))
        for member_index,(member,node) in enumerate(zip(record["q_star"],xml_stars),start=1):
            label = f"q{record['q_index']}/star{member_index}"
            raw = xml_force(node,nat)
            q_cart = xml_numbers(node.find("Q_POINT"))
            ledger.check(f"{label} independent XML force ordering", np.array_equal(raw,modes.force_tensor(member["force_matrix_ry_per_bohr2"])))
            ledger.check(f"{label} independent XML q vector", np.array_equal(q_cart,member["q_cartesian_2pi_alat"]))
            index = tuple((np.rint((at @ q_cart)*grid_shape).astype(int) % grid_shape).tolist())
            reference = replay_rows[index]
            q_error = np.abs(at @ reference["q_cartesian_2pi_alat"]-reference["q_fractional"])
            q_budget = (0.5*10**-6)*np.sum(np.abs(at),axis=1)+32*np.finfo(float).eps
            ledger.check(f"{label} native replay q-coordinate units", np.all(q_error <= q_budget))
            reported = np.sort(reference["frequency_cm_1"])*cm_1_to_thz
            reference_lambda = np.sign(reported)*(reported/modes.RY_TO_THZ)**2
            low,high = frequency_interval_qe(reported,half_width_thz)
            operator = modes.mass_weighted_operator(raw,masses,units="qe_ry")
            raw_residual = modes.operator_residuals(operator,masses)
            # Reproduce named native Hermitian processing, never ASR.
            hermitian = modes.project_hermitian(operator)
            spectrum = modes.signed_modes(hermitian,units="qe_ry")
            bound = spectral_budget(hermitian,reference_lambda,operation_scale)
            low,high = low-bound,high+bound
            frequency_bound = frequency_budget(spectrum.eigenvalues,reference_lambda,bound)+half_width_thz
            frequency_bound += 32*np.finfo(float).eps*np.abs(reported)
            absolute_errors = np.abs(spectrum.frequency_thz-reported)
            ledger.check(f"{label} serialized native signed squared spectrum",
                         np.all((spectrum.eigenvalues >= low) & (spectrum.eigenvalues <= high)))
            ledger.check(f"{label} absolute THz spectrum including acoustic residuals", np.all(absolute_errors <= frequency_bound))
            ledger.check(f"{label} eigensolver residual", spectrum.eigenpair_residual <= bound*np.sqrt(3*nat))
            ledger.check(f"{label} factor-two force-unit mutation exposed",
                         np.any((2*spectrum.eigenvalues < low) | (2*spectrum.eigenvalues > high)))
            rows.append({"q_index":record["q_index"],"star_index":member_index,
                         "q_cartesian_2pi_alat":member["q_cartesian_2pi_alat"],
                         "matrix_derived_signed_thz":spectrum.frequency_thz.tolist(),
                         "native_serialized_operator_signed_thz":reported.tolist(),
                         "ph_in_memory_signed_thz":ph_reported.tolist(),
                         "ph_in_memory_difference_thz":(spectrum.frequency_thz-ph_reported).tolist(),
                         "absolute_error_thz":absolute_errors.tolist(),
                         "absolute_budget_thz":frequency_bound.tolist(),
                         "absolute_eigenvalue_roundoff_budget_qe":bound,
                         "native_rounding_eigenvalue_interval_qe":[low.tolist(),high.tolist()],
                         "degenerate_index_clusters":degeneracy_clusters(reference_lambda,float(np.max(high-low))/2),
                         "raw_residuals_qe_units":raw_residual,
                         "explicit_hermitian_change_frobenius":float(np.linalg.norm(hermitian-operator)),
                         "eigenpair_residual_qe":spectrum.eigenpair_residual,
                         "orthogonality_residual":spectrum.orthogonality_residual})
    force_r = lattice.q_grid_to_force_lattice(grid)
    reconstruction = float(np.linalg.norm(lattice.force_lattice_to_q_grid(force_r)-grid))
    bound = 256*np.finfo(float).eps*np.linalg.norm(grid)
    ledger.check("raw pilot complete-grid transform reconstruction", reconstruction <= bound)
    return {"all_q_star_spectra":rows,"q_count":int(np.prod(grid.shape[:3])),
            "raw_lattice_residuals_ry_per_bohr2":lattice.lattice_residuals(force_r),
            "raw_grid_reconstruction_error":reconstruction,
            "native_frequency_scope":"Replay the serialized complete q grid through q2r and matdyn with no ASR/polar terms; compare every star member at its original sampled q.",
            "ph_frequency_scope":"Original PH OMEGA tags diagonalize its in-memory matrix, not the additional q2qstar written average; differences are retained as diagnostics, not identified with an eigensolver error.",
            "degeneracy_scope":"Sorted signed spectra and absolute intervals; no individual eigenvector comparison inside degenerate subspaces.",
            "physical_scope":"Unconverged Ni1 native engineering/unit fixture only. Gamma acoustic residuals remain raw, with no ASR, zero assignment, stability or material interpretation."}
