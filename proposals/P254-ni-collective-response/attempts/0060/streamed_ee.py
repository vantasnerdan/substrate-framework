"""Exact A46 EE, bounded in history length; q-additive numeric block ABI.

Original factors and an electron(a,b)->(greater,lesser) callback are required.
Each reverse leg is obtained by a separate callback, never inferred. Times are
exact values in callback index order (not assumed sorted). Static superoperators
are retained for a bounded external group and reused across every time tile.
No whole-history operands, Green bank, or nine-channel accumulator is allocated.
"""
from dataclasses import dataclass
from pathlib import Path
import hashlib
import json
import numpy as np
from all_slot_coulomb_shared import static_superoperators, exchange_superoperator
from all_slot_coulomb_blas import direct_q
from all_slot_coulomb_contractions import KEYS, validate_physical
STORED_KEYS = ('greater_direct','greater_exchange','lesser_direct','lesser_exchange')


@dataclass(frozen=True)
class Block:
    q: int
    external: tuple
    pairs: tuple
    values: dict


def _pair_tiles(nt, size, owned_pairs):
    # An iterator, not an nt^2 list. Explicit subsets preserve caller order.
    source = ((a,b) for a in range(nt) for b in range(nt)) if owned_pairs is None else iter(owned_pairs)
    tile = []
    for a,b in source:
        if not (0 <= a < nt and 0 <= b < nt):
            raise ValueError('time-pair index outside exact grid')
        tile.append((int(a),int(b)))
        if len(tile) == size:
            yield tuple(tile)
            tile = []
    if tile:
        yield tuple(tile)


def _operands(factors, q, electron, pairs):
    c = factors.c
    qm = int(c.minus[q])
    p = c.feature_partner[q][:len(factors.weights[q])]
    aa = c.target[qm,np.arange(c.ns)]
    raw = c.density[qm][:,p].conj()
    lq = c.density[q].conj()
    weight = (factors.weights[q]/c.n)[None,:,None,None]
    x = np.empty((2*len(pairs),c.ns,len(p),c.m,c.m),complex)
    y = np.empty_like(x)
    legs = []
    for k,(a,b) in enumerate(pairs):
        f,r = electron(a,b),electron(b,a)
        for value in (*f,*r):
            if value.shape != (c.ns,c.m,c.m):
                raise ValueError('electron callback must retain all explicit slots/bands')
        legs.append((f,r))
        for j,(ff,rr) in enumerate(((f[0],r[1]),(f[1],r[0]))):
            x[2*k+j] = (lq @ ff[:,None])*weight
            y[2*k+j] = rr[aa,None] @ raw @ ff[:,None]
    return x,y,legs


def stream_ee(factors, electron, times, consume, *, owned_q=None,
              owned_pairs=None, pair_batch=2, external_batch=2, static_slots=64,
              source_bindings=()):
    """Deliver exact q contributions, all 9 channels, shape [external,pair,m,m].

    consume is synchronous; it must persist/reduce blocks before returning and
    must not retain them to obtain bounded memory. q and ordered pair ownership
    are explicit: summing disjoint q partitions recovers the complete kernel.
    Within each q every requested pair and all external slots are emitted once.
    A subset is an explicitly partial kernel, never represented as full physics.
    Return manifest after successful synchronous delivery; exceptions propagate.
    No restart/prefix semantics are implied by this evaluator.
    """
    c = factors.c
    times = np.asarray(times,dtype=float)
    if times.ndim != 1 or not len(times) or not np.isfinite(times).all():
        raise ValueError('nonempty finite exact time grid required')
    if any(type(v) is not int or v <= 0 for v in (pair_batch,external_batch,static_slots)):
        raise ValueError('positive integer tile dimensions required')
    qs = tuple(range(c.nq)) if owned_q is None else tuple(owned_q)
    if len(set(qs)) != len(qs) or any(type(q) is not int or not 0 <= q < c.nq for q in qs):
        raise ValueError('unique owned outer q indices required')
    if owned_pairs is not None:
        owned_pairs = tuple(tuple(p) for p in owned_pairs)
        if len(set(owned_pairs)) != len(owned_pairs):
            raise ValueError('duplicate owned time pairs')
        # Validate before first numeric side effect, including empty subsets.
        for _ in _pair_tiles(len(times),pair_batch,owned_pairs):
            pass
    validate_physical(factors)
    manifest = dict(schema='p254-a60-streamed-ee-v1',
        times=times.tolist(), time_pair_order='a-major-b-minor; greater then lesser operands',
        owned_pairs='all' if owned_pairs is None else owned_pairs,
        owned_outer_q=qs, global_outer_q=c.nq, exchange_inner_h=c.nq,
        external_slots=c.ns, bands=c.m, channels=KEYS,
        layout='external,pair,band,band; complex128', additive_axis='outer_q',
        source_bindings=list(source_bindings), pair_batch=pair_batch,
        external_batch=external_batch, static_slots=static_slots,
        lifetime='consume synchronously; values valid for callback duration',
        reverse_legs='independent electron(b,a) callback',
        full_physical_scope=(len(qs)==c.nq and owned_pairs is None), blocks=0,
        static_payload_peak=0, operands_payload_peak=0, output_payload_peak=0)
    for q in qs:
        for begin in range(0,c.ns,static_slots):
            group = np.arange(begin,min(begin+static_slots,c.ns))
            matrix,b,cc = static_superoperators(factors,q,group)
            manifest['static_payload_peak'] = max(manifest['static_payload_peak'],matrix.nbytes)
            for pairs in _pair_tiles(len(times),pair_batch,owned_pairs):
                x,y,legs = _operands(factors,q,electron,pairs)
                manifest['operands_payload_peak'] = max(manifest['operands_payload_peak'],x.nbytes+y.nbytes)
                # Bubble/direct work shares all slots, then only selected outputs
                # survive. Repeated per bounded static group, never per external.
                direct = np.empty((2,len(group),len(pairs),c.m,c.m),complex)
                for k,(f,r) in enumerate(legs):
                    direct[0,:,k] = direct_q(factors,q,f[0],r[1])[group]
                    direct[1,:,k] = direct_q(factors,q,f[1],r[0])[group]
                for start in range(0,len(group),external_batch):
                    end = min(start+external_batch,len(group))
                    sl = slice(start*c.nq,end*c.nq)
                    xb = x[:,b[sl]].transpose(1,0,2,3,4)
                    yb = y[:,cc[sl]].transpose(1,0,2,3,4)
                    value = exchange_superoperator(matrix[sl],xb,yb)
                    value = value.reshape(end-start,c.nq,2*len(pairs),c.m,c.m).sum(axis=1)
                    values = {}
                    mask = np.asarray([times[a]>=times[z] for a,z in pairs])[None,:,None,None]
                    for j,order in enumerate(('greater','lesser')):
                        values[order+'_direct'] = direct[j,start:end].copy()
                        values[order+'_exchange'] = value[:,j::2].copy()
                        values[order+'_total'] = values[order+'_direct']+values[order+'_exchange']
                    for channel in ('direct','exchange','total'):
                        values['retarded_'+channel] = (values['greater_'+channel]-values['lesser_'+channel])*mask
                    manifest['output_payload_peak'] = max(manifest['output_payload_peak'],sum(v.nbytes for v in values.values()))
                    consume(Block(q,tuple(group[start:end].tolist()),pairs,values))
                    manifest['blocks'] += 1
                    del xb,yb,value,values
                del x,y,legs,direct
            del matrix,b,cc
    return manifest


class NumericSink:
    """Append-only numeric sink; no nt² array. close commits an offset manifest.

    Reader yields persisted blocks, retaining only one. Four independent channels
    are stored in STORED_KEYS order; total/retarded are reconstructed by the same
    arithmetic and exact saved time mask. Stream, index and manifest are required.
    """
    def __init__(self,path):
        self.path = Path(path)
        self.file = self.path.open('xb')
        self.index_path = Path(str(path)+'.index.jsonl')
        self.index = self.index_path.open('x')
        self.count = 0
    def __call__(self,block):
        offset = self.file.tell()
        np.save(self.file,np.stack([block.values[k] for k in STORED_KEYS]),allow_pickle=False)
        self.index.write(json.dumps(dict(offset=offset,end=self.file.tell(),q=block.q,
            external=block.external,pairs=block.pairs))+'\n')
        self.count += 1
    def close(self,manifest):
        self.file.close()
        self.index.close()
        h = hashlib.sha256()
        with self.path.open('rb') as f:
            for chunk in iter(lambda:f.read(1048576),b''):
                h.update(chunk)
        ih = hashlib.sha256()
        with self.index_path.open('rb') as f:
            for chunk in iter(lambda:f.read(1048576),b''):
                ih.update(chunk)
        record = dict(manifest=manifest,numeric=dict(path=str(self.path.resolve()),
            bytes=self.path.stat().st_size,sha256=h.hexdigest()),
            index=dict(path=str(self.index_path.resolve()),bytes=self.index_path.stat().st_size,
                sha256=ih.hexdigest()),block_count=self.count,stored_channels=STORED_KEYS)
        self.path.with_suffix(self.path.suffix+'.json').write_text(json.dumps(record,indent=2)+'\n')
        return record


def read_blocks(path):
    record = json.loads(Path(str(path)+'.json').read_bytes())
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):
            h.update(chunk)
        if h.hexdigest()!=record['numeric']['sha256'] or f.tell()!=record['numeric']['bytes']:
            raise ValueError('numeric stream binding differs')
        index = Path(record['index']['path'])
        ih = hashlib.sha256()
        with index.open('rb') as source:
            for chunk in iter(lambda:source.read(1048576),b''):
                ih.update(chunk)
            if ih.hexdigest()!=record['index']['sha256'] or source.tell()!=record['index']['bytes']:
                raise ValueError('block index binding differs')
        count = 0
        with index.open() as source:
            for line in source:
                item = json.loads(line)
                f.seek(item['offset'])
                array = np.load(f,allow_pickle=False)
                if f.tell()!=item['end'] or array.shape!=(len(STORED_KEYS),len(item['external']),len(item['pairs']),record['manifest']['bands'],record['manifest']['bands']):
                    raise ValueError('numeric block extent/layout differs')
                values = dict(zip(STORED_KEYS,array))
                times = record['manifest']['times']
                mask = np.asarray([times[a]>=times[b] for a,b in item['pairs']])[None,:,None,None]
                for order in ('greater','lesser'):
                    values[order+'_total'] = values[order+'_direct']+values[order+'_exchange']
                for channel in ('direct','exchange','total'):
                    values['retarded_'+channel] = (values['greater_'+channel]-values['lesser_'+channel])*mask
                count += 1
                yield Block(item['q'],tuple(item['external']),tuple(map(tuple,item['pairs'])),values)
        if count != record['block_count']:
            raise ValueError('block index count differs')
