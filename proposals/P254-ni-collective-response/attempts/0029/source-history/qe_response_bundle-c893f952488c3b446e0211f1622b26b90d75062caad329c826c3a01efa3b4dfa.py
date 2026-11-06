"""Verified QE response bundles; no ndarray is decoded before closure verification.

Native identity, dimension and index strings are parsed with Python's arbitrary
precision integers. Arrays keep their actual Fortran ranks, including FFT padding.
All electronic/PW and external-operator documents are retained for a future
full-complement, two-branch resolvent; this loader does not supply that driver.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

_DECIMAL = re.compile(r"-?(0|[1-9][0-9]*)\Z")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_DTYPES = {"float64": "f8", "complex128": "c16", "int32": "i4"}


def native_int(value: Any, label: str) -> int:
    """Require the lossless native decimal-string representation."""
    if not isinstance(value, str) or not _DECIMAL.fullmatch(value):
        raise ValueError(f"{label}: expected a native decimal integer string")
    return int(value)


def _unique_object(items: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _json(data: bytes) -> dict[str, Any]:
    value = json.loads(data, object_pairs_hook=_unique_object,
                       parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))
    if not isinstance(value, dict):
        raise ValueError("bundle metadata must be a JSON object")
    return value


def _basename(value: Any) -> str:
    if not isinstance(value, str) or value in ("", ".", "..") or Path(value).name != value:
        raise ValueError("bundle paths must be relative basenames")
    return value


def _metadata(value: Any) -> Any:
    if isinstance(value, str) and _DECIMAL.fullmatch(value):
        return int(value)
    if isinstance(value, dict):
        return MappingProxyType({key: _metadata(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_metadata(item) for item in value)
    return value


@dataclass(frozen=True)
class DatasetDescriptor:
    name: str
    path: str
    dtype: str
    shape: tuple[int, ...]
    units: str
    index_base: int | None
    bytes: int
    sha256: str


@dataclass(frozen=True)
class ResponseDocument:
    path: str
    metadata: Mapping[str, Any]
    descriptors: Mapping[str, DatasetDescriptor]
    datasets: Mapping[str, NDArray[Any]]

    def array(self, name: str, shape: tuple[int, ...] | None = None,
              dtype: str | None = None) -> NDArray[Any]:
        array = self.datasets[name]
        descriptor = self.descriptors[name]
        if shape is not None and array.shape != shape:
            raise ValueError(f"{self.path}:{name}: shape {array.shape}, expected {shape}")
        if dtype is not None and descriptor.dtype != dtype:
            raise ValueError(f"{self.path}:{name}: expected {dtype}, got {descriptor.dtype}")
        return array


@dataclass(frozen=True)
class QEResponseBundle:
    root: Path
    background: ResponseDocument
    states: tuple[ResponseDocument, ...]
    bare: tuple[ResponseDocument, ...]
    documents: Mapping[str, ResponseDocument]
    manifest: Mapping[str, Any]
    buffers: Mapping[str, bytes]


def load_bundle(path: str | Path) -> QEResponseBundle:
    """Load one per-q directory closed by ``manifest.json``.

    The manifest's files list binds every metadata JSON and binary payload with
    {path, bytes: decimal-string, sha256}; closed=true and byte_order are required.
    Metadata objects own a datasets mapping. Each descriptor declares native
    dtype, decimal-string shape, order='F', units and index_base='1' for indices.
    Descriptor bytes/SHA may also be present, and must agree with the closure.
    Verification uses the exact same bytes subsequently decoded (no reread race).
    """
    root = Path(path).resolve()
    manifest_bytes = (root / "manifest.json").read_bytes()
    manifest = _json(manifest_bytes)
    if manifest.get("schema_version") != 1 or manifest.get("closed") is not True:
        raise ValueError("expected a closed version-1 response manifest")
    if manifest.get("byte_order") not in ("little", "big"):
        raise ValueError("manifest must declare actual native byte_order")
    records = manifest.get("files")
    if not isinstance(records, list) or not records:
        raise ValueError("manifest must bind its complete files list")
    raw: dict[str, bytes] = {}
    bindings: dict[str, tuple[int, str]] = {}
    # Complete closure verification precedes dtype interpretation or reshape.
    for record in records:
        name = _basename(record["path"])
        if name == "manifest.json" or name in raw:
            raise ValueError(f"duplicate or self-referential manifest file: {name}")
        size = native_int(record["bytes"], f"{name}.bytes")
        digest = record["sha256"]
        if size < 0 or not isinstance(digest, str) or not _SHA256.fullmatch(digest):
            raise ValueError(f"invalid content binding: {name}")
        data = (root / name).read_bytes()
        if len(data) != size or hashlib.sha256(data).hexdigest() != digest:
            raise ValueError(f"closed payload byte count/SHA mismatch: {name}")
        raw[name] = data
        bindings[name] = (size, digest)
    documents: dict[str, ResponseDocument] = {}
    referenced: set[str] = set()
    endian = "<" if manifest["byte_order"] == "little" else ">"
    for name, data in raw.items():
        if not name.endswith(".json"):
            continue
        meta = _json(data)
        records_by_name = meta.get("datasets")
        if not isinstance(records_by_name, dict):
            raise ValueError(f"{name}: datasets mapping required")
        descriptors: dict[str, DatasetDescriptor] = {}
        arrays: dict[str, NDArray[Any]] = {}
        for semantic, record in records_by_name.items():
            payload = _basename(record["path"])
            if payload not in raw or payload.endswith(".json"):
                raise ValueError(f"{name}:{semantic}: unbound binary payload")
            dtype = record["dtype"]
            if dtype not in _DTYPES or record.get("order") != "F":
                raise ValueError(f"{name}:{semantic}: unsupported native dtype/order")
            shape_value = record["shape"]
            if not isinstance(shape_value, list) or not shape_value:
                raise ValueError(f"{name}:{semantic}: actual native rank required")
            shape = tuple(native_int(n, f"{semantic}.shape") for n in shape_value)
            if any(n < 0 for n in shape):
                raise ValueError(f"{name}:{semantic}: negative dimension")
            units = record["units"]
            if not isinstance(units, str) or not units:
                raise ValueError(f"{name}:{semantic}: units required")
            index_base = (native_int(record["index_base"], f"{semantic}.index_base")
                          if "index_base" in record else None)
            if index_base is not None and index_base != 1:
                raise ValueError("only native one-based indices are supported")
            size, digest = bindings[payload]
            for key, expected in (("bytes", size), ("sha256", digest)):
                if key in record:
                    actual = native_int(record[key], f"{semantic}.bytes") if key == "bytes" else record[key]
                    if actual != expected:
                        raise ValueError(f"{name}:{semantic}: descriptor disagrees with closure")
            np_dtype = np.dtype(endian + _DTYPES[dtype])
            if math.prod(shape) * np_dtype.itemsize != size:
                raise ValueError(f"{name}:{semantic}: dtype/shape byte count mismatch")
            array = np.frombuffer(raw[payload], dtype=np_dtype).reshape(shape, order="F")
            if dtype != "int32" and not np.isfinite(array).all():
                raise ValueError(f"{name}:{semantic}: nonfinite native data")
            arrays[semantic] = array
            descriptors[semantic] = DatasetDescriptor(semantic, payload, dtype, shape,
                                                       units, index_base, size, digest)
            referenced.add(payload)
        documents[name] = ResponseDocument(name, _metadata({k: v for k, v in meta.items()
                                                            if k != "datasets"}),
                                            MappingProxyType(descriptors), MappingProxyType(arrays))
    if "background.json" not in documents:
        raise ValueError("closed bundle missing background.json")
    if set(raw) != set(documents) | referenced:
        raise ValueError("closure contains unowned/non-dataset binary payloads")
    unknown = [name for name in documents if name != "background.json" and
               not name.startswith(("state-", "bare-"))]
    if unknown:
        raise ValueError(f"unknown response document roles: {unknown}")
    bundle = QEResponseBundle(root, documents["background.json"],
                              tuple(doc for name, doc in sorted(documents.items()) if name.startswith("state-")),
                              tuple(doc for name, doc in sorted(documents.items()) if name.startswith("bare-")),
                              MappingProxyType(documents), _metadata(manifest),
                              MappingProxyType({"manifest.json": manifest_bytes, **raw}))
    from .qe_response_schema import validate_bundle
    validate_bundle(bundle)
    return bundle
