"""Actual symmetric-preparation compression, never a principal raw-D block."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class ActiveBasis:
    columns: np.ndarray
    rho0: np.ndarray
    labels: np.ndarray
    selected: np.ndarray
    diagnostics: dict


def prepare_active_basis(psi, raw_d, labels, interval, bands=12, *, gram=None, populations=None):
    psi, raw_d, labels = np.asarray(psi), np.asarray(raw_d), np.asarray(labels)
    if psi.ndim != 2 or raw_d.shape != (psi.shape[1],)*2 or labels.shape != (psi.shape[1],):
        raise ValueError('full original raw frame, D and saved energy labels required')
    if bands not in (8, 12) or bands > len(labels):
        raise ValueError('nested saved-band size must be 8 or 12')
    s = psi.conj().T @ psi if gram is None else np.asarray(gram)
    if s.shape != raw_d.shape or not all(np.isfinite(a).all() for a in (s, raw_d, labels)):
        raise ValueError('finite full original Gram/D/labels required')
    values, vectors = np.linalg.eigh(s)
    if np.min(values) <= 0:
        raise ValueError('original Gram is not positive definite; no repair permitted')
    root = (vectors*np.sqrt(values)) @ vectors.conj().T
    inverse = (vectors/np.sqrt(values)) @ vectors.conj().T
    rho = root @ raw_d @ root
    low, high = map(float, interval)
    if not np.isfinite([low, high]).all() or low > high:
        raise ValueError('saved neutral chemical-potential interval required')
    distance = np.maximum(np.maximum(low-labels, labels-high), 0.)
    order = np.lexsort((np.arange(len(labels)), labels, distance))
    selected = order[:bands]
    remainder = order[bands:]
    q = psi @ inverse[:, selected]
    active = rho[np.ix_(selected, selected)]
    offdiag = rho-np.diag(np.diag(rho))
    diagnostics = dict(
        raw_gram_hermiticity=float(np.linalg.norm(s-s.conj().T)),
        raw_d_hermiticity=float(np.linalg.norm(raw_d-raw_d.conj().T)),
        gram_minimum=float(values.min()), gram_condition=float(values.max()/values.min()),
        active_orthonormality=float(np.linalg.norm(q.conj().T @ q-np.eye(bands))),
        preparation_offdiagonal=float(np.linalg.norm(offdiag)),
        active_passive_coherence=float(np.linalg.norm(rho[np.ix_(selected, remainder)])),
        original_population_trace_real=float(np.trace(rho).real),
        original_population_trace_imag=float(np.trace(rho).imag),
        active_population_trace_real=float(np.trace(active).real),
        active_population_trace_imag=float(np.trace(active).imag),
        raw_to_prepared_coordinate_residual=float(np.linalg.norm(root @ inverse-np.eye(len(labels)))),
        selected_labels_ry=labels[selected].tolist(), selected_raw_labels=selected.tolist(),
        neutral_interval_ry=[low, high], no_population_refill=True)
    if populations is not None:
        populations = np.asarray(populations)
        if populations.shape != (len(labels),):
            raise ValueError('full saved population ledger required')
        diagnostics['population_coordinate_defect'] = float(np.linalg.norm(rho-np.diag(populations)))
    # The settled preparation is population-diagonal in this frame. A material
    # violation is rejected, not hidden by choosing a principal block or refilling.
    tolerance = 1e-8*max(1., float(np.linalg.norm(rho)))
    if (diagnostics['preparation_offdiagonal'] > tolerance
            or diagnostics.get('population_coordinate_defect', 0.) > tolerance
            or diagnostics['active_orthonormality'] > 1e-8*bands):
        raise ValueError('saved population-frame identity fails: '+repr(diagnostics))
    return ActiveBasis(q, active, labels[selected], selected, diagnostics)
