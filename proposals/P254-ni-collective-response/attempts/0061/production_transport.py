"""A61 result-only transport; unchanged A54 physics, caller-owned mapped banks.

No NPZ/auxiliary decoded banks. A slot commits after numeric fsync+digest; only
its declared newly-created packed staging may then be unlinked. Resume is at a
committed slot, never inside the fixed-history native solve. Diagnostic consumers
continue to use the immutable A54 supplied_keldysh.py companion interface.
"""
from pathlib import Path
import hashlib
import json
import os
import struct
import io
from numbers import Integral
import sys
import numpy as np

A49 = Path(__file__).resolve().parent.parent / '0049'
if str(A49) not in sys.path:
    sys.path.append(str(A49))
from command_receipts import binding, run_recorded

HEADER = struct.Struct('<8siiiid')
SCHEMA = 'p254-a61-result-bank-v1'
SOURCEKIND = 'single-physical-Hermitian-contour'
def _byte_limit(value, size, label):
    if value is not None:
        if isinstance(value, bool) or not isinstance(value, Integral) or value <= 0:
            raise ValueError(label+' limit must be a positive integer or None')
        if size > value:
            raise ValueError(label+' exceeds caller-declared byte limit')


def plan_layout(time_count, dimension, slots=1, *, max_input_bytes=None,
                max_result_bytes=None, max_bank_bytes=None):
    """Exact shape/byte planning, no physical arrays or files allocated.

    Optional limits belong to the caller. None imposes only native/NumPy integer
    representability, never this campaign's 64MiB benchmark allocation.
    """
    values = (time_count, dimension, slots)
    if any(isinstance(v, bool) or not isinstance(v, Integral) for v in values):
        raise ValueError('layout sizes must be integers')
    count, d, slots = map(int, values)
    if count < 2 or count-1 > np.iinfo(np.int32).max or d < 1 or d*d > np.iinfo(np.int32).max or slots < 1:
        raise ValueError('layout outside native header/matrix integer representation')
    cells = count*(count+1)//2
    shape = (slots,count,count,d,d)
    array_bytes = slots*count*count*d*d*16
    header = io.BytesIO()
    np.lib.format.write_array_header_1_0(header,dict(descr='<c16',fortran_order=False,shape=shape))
    input_bytes = HEADER.size+count*8+cells*d*d*64
    result_bytes = HEADER.size+count*8+cells*d*d*32
    bank_bytes = 2*(array_bytes+header.tell())
    if max(input_bytes,result_bytes,bank_bytes) > np.iinfo(np.intp).max:
        raise ValueError('layout exceeds addressable file/mapping representation')
    _byte_limit(max_input_bytes,input_bytes,'input transport')
    _byte_limit(max_result_bytes,result_bytes,'result transport')
    _byte_limit(max_bank_bytes,bank_bytes,'final numeric bank')
    return dict(shape=list(shape),triangle_cells=cells,input_bytes=input_bytes,
                result_bytes=result_bytes,array_payload_bytes=array_bytes,
                npy_header_bytes=header.tell(),bank_bytes=bank_bytes)


def atomic_json(path, value):
    path = Path(path)
    temp = path.with_name(path.name + '.new')
    with temp.open('w') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp, path)
    fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def array_digest(array):
    """Canonical C-order little-endian numeric hash, at most one dxd copy."""
    digest = hashlib.sha256()
    for index in np.ndindex(array.shape[:-2]):
        digest.update(np.asarray(array[index], dtype='<c16', order='C').tobytes())
    return dict(shape=list(array.shape), dtype='<c16', sha256=digest.hexdigest(), bytes=array.size*16)


def validate(times, arrays, statistics, order, beta, *, max_input_bytes=None, max_result_bytes=None):
    raw_times = np.asarray(times)
    if raw_times.dtype.kind != 'f' or raw_times.dtype.itemsize != 8:
        raise ValueError('times require float64; no silent precision conversion')
    times = np.asarray(raw_times, dtype='<f8')
    if (times.ndim != 1 or len(times) <= order or order not in range(1, 6)
            or statistics not in (-1, 1) or not np.isfinite(beta) or beta <= 0):
        raise ValueError('invalid time/order/statistics/beta')
    dt = times[1]-times[0]
    if times[0] != 0 or dt <= 0 or not np.isfinite(times).all() or not np.allclose(times, np.arange(len(times))*dt, rtol=1e-12, atol=1e-12):
        raise ValueError('uniform finite time history from zero required; no interpolation')
    shape = arrays[0].shape
    if len(shape) != 4 or shape[:2] != (len(times),len(times)) or shape[2] < 1 or shape[2] != shape[3]:
        raise ValueError('expected [time,time,row,column]')
    plan_layout(len(times),shape[2],max_input_bytes=max_input_bytes,max_result_bytes=max_result_bytes)
    for index, a in enumerate(arrays):
        if a.shape != shape or a.dtype.kind != 'c' or a.dtype.itemsize != 16:
            raise ValueError('all components require equal shapes and complex128')
        scale, defect = 0., 0.
        for n in range(len(times)):
            for m in range(len(times)):
                block = a[n,m]
                if not np.isfinite(block).all():
                    raise ValueError('nonfinite component')
                if index % 2:
                    scale = max(scale, float(np.max(np.abs(block))))
                    defect = max(defect, float(np.max(np.abs(block+a[m,n].conj().T))))
                elif m > n and np.any(block != 0):
                    raise ValueError('retarded upper triangle must be exactly zero')
        if index % 2 and defect > 1e-10*(1+scale):
            raise ValueError('independent lesser legs unsupported; no projection/repair')
    return times, shape[2]


class ResultView:
    """Read-only zero-copy packed triangles, valid only until close()."""
    def __init__(self, path, header=None, times=None, *, max_result_bytes=None):
        self.path = Path(path)
        with self.path.open('rb') as stream:
            raw = stream.read(HEADER.size)
        if len(raw) != HEADER.size or (header is not None and raw != header):
            raise ValueError('result header mismatch')
        magic, nt, d, order, statistics, beta = HEADER.unpack(raw)
        if magic != b'P254VIE1' or nt < order or order not in range(1,6) or d < 1 or statistics not in (-1,1) or not np.isfinite(beta) or beta <= 0:
            raise ValueError('invalid minimal result header')
        count, cells = nt+1, (nt+1)*(nt+2)//2
        expected = plan_layout(count,d,max_result_bytes=max_result_bytes)['result_bytes']
        if self.path.stat().st_size != expected:
            raise ValueError('minimal result extent mismatch; diagnostic six-component files unsupported')
        self._mapped = np.memmap(self.path, mode='r', dtype='u1')
        self.times = np.ndarray((count,), dtype='<f8', buffer=self._mapped, offset=HEADER.size)
        if (times is not None and not np.array_equal(self.times,times)) or not np.isfinite(self.times).all() or self.times[0] != 0 or self.times[1] <= 0 or not np.allclose(self.times,np.arange(count)*self.times[1],rtol=1e-12,atol=1e-12):
            self.close()
            raise ValueError('result time mismatch/nonuniform times')
        self.ret = np.ndarray((cells,d,d),dtype='<c16',buffer=self._mapped,offset=HEADER.size+count*8)
        self.less = np.ndarray((cells,d,d),dtype='<c16',buffer=self._mapped,offset=HEADER.size+count*8+cells*d*d*16)
        self.statistics, self.order, self.beta = statistics, order, beta

    def copy_into(self, ret, less):
        count, d = len(self.times), self.ret.shape[-1]
        if ret.shape != (count,count,d,d) or less.shape != ret.shape or ret.dtype != np.dtype('<c16') or less.dtype != ret.dtype:
            raise ValueError('destination shape/dtype mismatch')
        # One physical Hermitian contour ABI, not manufacture of paired-q data.
        offset = 0
        for n in range(count):
            for m in range(count):
                ret[n,m] = 0
            for m in range(n+1):
                r, l = self.ret[offset+m], self.less[offset+m]
                if not np.isfinite(r).all() or not np.isfinite(l).all():
                    raise ValueError('nonfinite result')
                ret[n,m], less[m,n] = r, l
                if m != n:
                    less[n,m] = -l.conj().T
                elif np.max(np.abs(l+l.conj().T)) > 1e-10*(1+np.max(np.abs(l))):
                    raise ValueError('raw result lesser diagonal exchange failure')
            offset += n+1

    def close(self):
        self.ret = self.less = self.times = None
        if getattr(self, '_mapped', None) is not None:
            self._mapped._mmap.close()
            self._mapped = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


class FinalBank:
    """Single-writer final numeric bank with durable per-slot recovery receipts.

    Call finalize() before downstream full-bank use. load_solution() requires the
    final receipt. Resume validates all committed slot bytes before returning.
    No retained array view may be used after close(); context lifetime is explicit.
    """
    def __init__(self, directory, times=None, dimension=None, slots=None, *, statistics=None, order=None, beta=1., resume=False, max_bank_bytes=None):
        self.path = Path(directory)
        if resume:
            self.meta = json.loads((self.path/'bank.json').read_text())
            if self.meta['schema'] != SCHEMA:
                raise ValueError('unsupported bank schema')
        else:
            times = np.asarray(times)
            if (times.dtype != np.dtype('<f8') or times.ndim != 1 or len(times) < 2 or dimension is None or dimension < 1 or slots is None or slots < 1 or order not in range(1,6) or len(times) <= order or statistics not in (-1,1) or not np.isfinite(beta) or beta <= 0):
                raise ValueError('invalid bank layout')
            if times[0] != 0 or times[1] <= 0 or not np.isfinite(times).all() or not np.allclose(times,np.arange(len(times))*times[1],rtol=1e-12,atol=1e-12):
                raise ValueError('invalid bank times')
            layout = plan_layout(len(times),dimension,slots,max_bank_bytes=max_bank_bytes)
            shape = tuple(layout['shape'])
            self.path.mkdir(parents=True,exist_ok=False)
            (self.path/'receipts').mkdir()
            (self.path/'.staging').mkdir()
            self.meta = dict(schema=SCHEMA,sourcekind=SOURCEKIND,shape=list(shape),times=times.tolist(),statistics=statistics,order=order,beta=beta,layout=layout,max_bank_bytes=max_bank_bytes)
            atomic_json(self.path/'bank.json',self.meta)
            for name in ('ret','less'):
                a = np.lib.format.open_memmap(self.path/(name+'.npy'),mode='w+',dtype='<c16',shape=shape)
                a.flush()
                a._mmap.close()
        self.ret = np.load(self.path/'ret.npy',mmap_mode='r+',allow_pickle=False)
        self.less = np.load(self.path/'less.npy',mmap_mode='r+',allow_pickle=False)
        if list(self.ret.shape) != self.meta['shape'] or self.less.shape != self.ret.shape or self.ret.dtype != np.dtype('<c16') or self.less.dtype != self.ret.dtype:
            raise ValueError('bank numeric layout mismatch')
        self.layout = plan_layout(self.ret.shape[1],self.ret.shape[-1],self.ret.shape[0],max_bank_bytes=max_bank_bytes)
        self.times = np.asarray(self.meta['times'],dtype='<f8')
        self.completed = {}
        for p in sorted((self.path/'receipts').glob('slot-*.json')):
            r = json.loads(p.read_text())
            slot = r['slot']
            if r['ret'] != array_digest(self.ret[slot]) or r['less'] != array_digest(self.less[slot]):
                raise ValueError('committed slot numeric digest mismatch')
            self.completed[slot] = r

    def commit(self, slot, provenance):
        for name, a in (('ret',self.ret),('less',self.less)):
            a.flush()
            with (self.path/(name+'.npy')).open('rb') as stream:
                os.fsync(stream.fileno())
        receipt = dict(slot=slot,ret=array_digest(self.ret[slot]),less=array_digest(self.less[slot]),provenance=provenance)
        atomic_json(self.path/'receipts'/('slot-%03d.json'%slot),receipt)
        self.completed[slot] = receipt
        return receipt

    def finalize(self):
        if set(self.completed) != set(range(self.ret.shape[0])):
            raise ValueError('bank has uncommitted slots')
        receipt = dict(schema=SCHEMA,bank=binding(self.path/'bank.json'),ret=binding(self.path/'ret.npy'),less=binding(self.path/'less.npy'),slots=[binding(self.path/'receipts'/('slot-%03d.json'%slot)) for slot in sorted(self.completed)])
        atomic_json(self.path/'complete.json',receipt)
        return self.path

    def close(self):
        for name in ('ret','less'):
            a = getattr(self,name,None)
            if a is not None:
                a._mmap.close()
                setattr(self,name,None)

    def __enter__(self):
        return self

    def __exit__(self,*exc):
        self.close()


def solve_supplied(times, source_ret, source_less, kernel_ret, kernel_less, *, statistics, order, binary, bank, slot, source_bindings, beta=1., timeout=120, max_input_bytes=None, max_result_bytes=None):
    """Stream one channel/slot into bank; return durable slot receipt, not arrays.

    source_bindings names actual original material/source files. Packed bytes and
    canonical full-array numeric hashes also bind this exact invocation. Failure
    retains declared staging and receipts; it never deletes a pre-existing path.
    """
    if not source_bindings:
        raise ValueError('actual source bindings required')
    if not 0 <= slot < bank.ret.shape[0]:
        raise ValueError('slot outside bank')
    if slot in bank.completed or (bank.path/'complete.json').exists():
        raise ValueError('slot/bank already committed; no overwrite')
    arrays = [np.asarray(a) for a in (source_ret,source_less,kernel_ret,kernel_less)]
    times,d = validate(times,arrays,statistics,order,beta,max_input_bytes=max_input_bytes,max_result_bytes=max_result_bytes)
    if (not np.array_equal(times,bank.times) or d != bank.ret.shape[-1] or statistics != bank.meta['statistics'] or order != bank.meta['order'] or beta != bank.meta['beta']):
        raise ValueError('solver/bank layout mismatch')
    owned = bank.path/'.staging'/('slot-%03d'%slot)
    owned.mkdir(exist_ok=False)
    packed, output = owned/'input.bin', owned/'result.bin'
    header = HEADER.pack(b'P254VIE1',len(times)-1,d,order,statistics,beta)
    with packed.open('xb') as stream:
        stream.write(header)
        stream.write(times.tobytes())
        for index,a in enumerate(arrays):
            for n in range(len(times)):
                block = a[:n+1,n] if index%2 else a[n,:n+1]
                stream.write(np.asarray(block,dtype='<c16',order='C').tobytes())
        stream.flush()
        os.fsync(stream.fileno())
    provenance = dict(algorithm='p254-a54-two-sided-keldysh-v1-unchanged',sourcekind=SOURCEKIND,sources=source_bindings,components=[array_digest(a) for a in arrays],input=binding(packed),binary=binding(binary),consumer=binding(__file__),staging_declared=[str(packed),str(output)],recovery='committed-slot-only; failed/in-flight staging retained; no silent replay')
    provenance['layout'] = plan_layout(len(times),d)
    provenance['caller_limits'] = dict(max_input_bytes=max_input_bytes,max_result_bytes=max_result_bytes,max_bank_bytes=bank.meta.get('max_bank_bytes'))
    atomic_json(owned/'invocation.json',provenance)
    argv = [str(Path(binary).resolve()),str(packed.resolve()),str(output.resolve())]
    if max_input_bytes is not None:
        argv.append(str(max_input_bytes))
    run_recorded(['/usr/bin/time','-v','-o',str((owned/'native.time').resolve()),*argv],owned,'native',timeout)
    provenance.update(output=binding(output),native=binding(owned/'native.json'))
    with ResultView(output,header,times,max_result_bytes=max_result_bytes) as result:
        result.copy_into(bank.ret[slot],bank.less[slot])
    receipt = bank.commit(slot,provenance)
    # Receipt authenticates exact final bytes before releasing only owned staging.
    for p in (packed,output):
        p.unlink()
    atomic_json(owned/'released.json',dict(committed=binding(bank.path/'receipts'/('slot-%03d.json'%slot)),released=provenance['staging_declared']))
    return receipt


def load_solution(path):
    """Return read-only mapped final-bank components; caller must close mappings."""
    path = Path(path)
    receipt = json.loads((path/'complete.json').read_text())
    if receipt['schema'] != SCHEMA:
        raise ValueError('unsupported bank completion schema')
    for name, filename in (('bank','bank.json'),('ret','ret.npy'),('less','less.npy')):
        if binding(path/filename) != receipt[name]:
            raise ValueError('final bank digest mismatch')
    for slot_receipt in receipt['slots']:
        if binding(path/'receipts'/Path(slot_receipt['path']).name) != slot_receipt:
            raise ValueError('final bank slot receipt mismatch')
    meta = json.loads((path/'bank.json').read_text())
    return dict(times=np.asarray(meta['times'],dtype='<f8'),statistics=meta['statistics'],order=meta['order'],beta=meta['beta'],ret=np.load(path/'ret.npy',mmap_mode='r',allow_pickle=False),less=np.load(path/'less.npy',mmap_mode='r',allow_pickle=False))
