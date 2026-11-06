"""Original physical density and shared bounded-PW projector moments.

The existing bridge owns density arithmetic. Cancellation-sensitive material
gradients require that canonical result; FFT sharing is not worth a second
density implementation. Original rectangular coefficients are still absorbed
on the smaller leg, with request-local immutable factors shared by the ionic
and kinetic contractions. Each column slice retires its factor/overlap caches.
"""
import numpy as np
from mixed_dyad_ionic_moments import ionic_gradient_block


def _factor_requests(terms):
    requests, work = {}, []
    for term in terms:
        nl, nr = term.left.shape[1], term.right.shape[1]
        width = min(nl, nr)
        if nr <= nl:
            left_key = ('left', id(term.left), id(term.D))
            right_key = ('plain', id(term.right))
            factors = ((left_key, term.left, term.D, 'left'),
                       (right_key, term.right, None, 'plain'))
        else:
            left_key = ('plain', id(term.left))
            right_key = ('right', id(term.right), id(term.D))
            factors = ((left_key, term.left, None, 'plain'),
                       (right_key, term.right, term.D, 'right'))
        for key, frame, coefficient, side in factors:
            if key not in requests or requests[key][3] < width:
                requests[key] = frame, coefficient, side, width
        work.append((term, left_key, right_key, width))
    return requests, work


def _column_blocks(requests, start, size):
    blocks = {}
    for key, (frame, coefficient, side, columns) in requests.items():
        if start >= columns:
            continue
        stop = min(start+size, columns)
        if side == 'plain':
            blocks[key] = frame[:, start:stop]
        elif side == 'left':
            blocks[key] = frame@coefficient[:, start:stop]
        else:
            blocks[key] = frame@coefficient[start:stop, :].conj().T
    return blocks




def _finish_gradient(coupling, channel, gradient):
    sampled = coupling.bridge.density_fft.to_g(channel)
    ref = coupling.inputs.reference
    for atom, name in enumerate(ref.atom_names):
        weighted = (sampled*coupling.vloc[name].conj()
                    *np.exp(1j*(coupling.gq@ref.positions_bohr[:, atom])))
        gradient[atom] += 1j*coupling.bridge.source_fft.omega*(weighted@coupling.gq)


def grid_moments(coupling, terms, *, field_momenta=None):
    """Return original grouped density/gradient and per-term kinetic L.

    Reuse the existing physical density solution, without reassociating its
    column reduction. Projector moments share bounded coefficient-absorbed PW
    factors while preserving each original dyad's finish and weighting order.
    """
    terms = tuple(terms)
    bridge, size = coupling.bridge, coupling.bridge.block_size
    requests, work = _factor_requests(terms)
    gradients = [np.zeros((coupling.nat, 3), complex) for _ in terms]
    diagonals = ([np.zeros(len(field_momenta), complex) for _ in terms]
                 if field_momenta is not None else None)
    largest = max((width for _, _, _, width in work), default=0)
    for start in range(0, largest, size):
        blocks = _column_blocks(requests, start, size)
        overlaps = {}
        for index, (term, left_key, right_key, width) in enumerate(work):
            if start >= width:
                continue
            count = min(size, width-start)
            left, right = blocks[left_key][:, :count], blocks[right_key][:, :count]
            gradients[index] += ionic_gradient_block(
                coupling, left, right, left_key, right_key, overlaps)
            if diagonals is not None:
                diagonals[index] += np.einsum('pi,pi->p', right.conj(), left)
    result, kinetic = {}, []
    for index, term in enumerate(terms):
        density = bridge.density_from_dyads(
            term.left, term.right, term.D, weight=term.weight)
        gradient = gradients[index]
        gradient *= term.weight
        _finish_gradient(coupling, density[:, bridge.ispin-1], gradient)
        if term.group not in result:
            result[term.group] = dict(
                density=np.zeros((bridge.source_fft.nnr, 2), complex),
                gradient=np.zeros(coupling.ndof, complex))
        result[term.group]['density'] += density
        result[term.group]['gradient'] += gradient.reshape(-1)
        kinetic.append(2*(field_momenta.T@diagonals[index])
                       if diagonals is not None else np.zeros(3, complex))
    return result, kinetic
