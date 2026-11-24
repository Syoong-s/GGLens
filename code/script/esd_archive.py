"""Versioned binary fit archives. The current fit points to one committed pair of files."""

from contextlib import contextmanager
from dataclasses import dataclass
import fcntl
from functools import cached_property
import hashlib
import math
import os
from pathlib import Path
import struct
from tempfile import mkdtemp

import numpy as np


VERSION = 1
CHAIN_MAGIC = b"FQCHAIN1"
CURVE_MAGIC = b"FQCURVE1"
CHAIN_HEADER = struct.Struct("<8sI32sQ4IQ4I")
CURVE_HEADER = struct.Struct("<8sI32s5I")
PARAMETERS = ("log10_M180m", "concentration", "f_cen", "R_sig_over_r180m", "b_h")


def read_exact(stream, count):
    if count > os.fstat(stream.fileno()).st_size - stream.tell():
        raise ValueError(f"Truncated fit archive: {stream.name}")
    value = stream.read(count)
    if len(value) != count:
        raise ValueError(f"Truncated fit archive: {stream.name}")
    return value


def pad_numeric_block(stream, writing=False):
    size = -stream.tell() % 8
    if writing:
        stream.write(bytes(size))
    elif read_exact(stream, size) != bytes(size):
        raise ValueError(f"Invalid fit archive padding: {stream.name}")


def map_array(stream, dtype, shape):
    dtype = np.dtype(dtype)
    offset = stream.tell()
    size = math.prod(shape) * dtype.itemsize
    end = offset + size
    if end > os.fstat(stream.fileno()).st_size:
        raise ValueError(f"Truncated fit archive array: {stream.name}")
    if size == 0:
        return np.empty(shape, dtype=dtype)
    result = np.memmap(stream, dtype=dtype, shape=shape, mode="r", offset=offset)
    stream.seek(end)
    return result


def expect_end(stream):
    if stream.read(1):
        raise ValueError(f"Unexpected trailing data in fit archive: {stream.name}")


def encode_chain(chain):
    """Keep the starting step of each state per walker, preserving all six double bit patterns."""
    tracks = chain.transpose(0, 2, 1, 3)
    bits = tracks.view(np.uint64)
    changed = np.empty(tracks.shape[:-1], dtype=bool)
    changed[..., 0] = True
    changed[..., 1:] = np.any(bits[..., 1:, :] != bits[..., :-1, :], axis=-1)
    states = np.empty((6, np.count_nonzero(changed)), dtype="<f8")
    for column in range(6):
        states[column] = tracks[..., column][changed]
    return np.packbits(changed, axis=-1, bitorder="little"), states


@dataclass
class ChainArchive:
    fit_id: bytes
    measurement_id: int
    measurement_path: str
    title: str
    shape: tuple
    burn_in: int
    seed: int
    mask: np.ndarray
    diagnostics: np.ndarray
    states: np.ndarray
    legacy_cosmology: np.ndarray

    @cached_property
    def starts(self):
        steps = self.shape[1]
        changed = np.unpackbits(self.mask, axis=-1, count=steps, bitorder="little")
        starts = np.flatnonzero(changed)
        if not np.all(changed[..., 0]) or len(starts) != self.states.shape[-1]:
            raise ValueError("Invalid MCMC state change bitmap")
        return starts

    @cached_property
    def weights(self):
        return np.diff(self.starts, append=math.prod(self.shape))

    def best(self):
        _, steps, walkers = self.shape
        likelihood = self.states[-1]
        candidates = np.flatnonzero(likelihood == np.max(likelihood))
        track, step = np.divmod(self.starts[candidates], steps)
        ensemble, walker = np.divmod(track, walkers)
        flat = (ensemble * steps + step) * walkers + walker
        winner = int(np.argmin(flat))
        return int(flat[winner]), float(likelihood[candidates[winner]])

    def expand(self):
        """Restore the original (ensemble, step, walker, parameter or log probability) array."""
        ensembles, steps, walkers = self.shape
        return np.repeat(self.states, self.weights, axis=1).reshape(
            6, ensembles, walkers, steps
        ).transpose(1, 3, 2, 0)


@dataclass
class CurveArchive:
    radius_count: int
    extra_radii: np.ndarray
    chain_indices: np.ndarray
    rows: np.ndarray
    best_row: int
    components: np.ndarray

    def radius(self, observed):
        result = np.union1d(self.extra_radii, observed)
        if len(result) != self.radius_count:
            raise ValueError("Fitted measurement does not match the profile radius grid")
        return result


def archive_id(chain, curves):
    signature = hashlib.sha256(struct.pack("<Q5I", chain.measurement_id, *chain.shape, chain.burn_in, chain.seed))
    title = chain.title.encode("utf-8")
    signature.update(struct.pack("<8I", chain.states.shape[-1], len(chain.legacy_cosmology), len(title),
                                 curves.radius_count, curves.best_row, len(curves.extra_radii),
                                 len(curves.chain_indices), curves.components.shape[-1]))
    signature.update(title)
    for values in (chain.mask, chain.legacy_cosmology, chain.diagnostics, chain.states,
                   curves.extra_radii, curves.chain_indices, curves.rows, curves.components):
        signature.update(memoryview(values).cast("B"))
    return signature.digest()


def read_chain(path):
    with path.open("rb") as stream:
        (magic, version, fit_id, source_id, ensembles, steps, walkers, burn, unique,
         seed, source_size, title_size, extra_cosmology) = CHAIN_HEADER.unpack(read_exact(stream, CHAIN_HEADER.size))
        samples = ensembles * steps * walkers
        if (magic != CHAIN_MAGIC or version != VERSION or min(ensembles, steps, walkers) < 1
                or samples > (2**31 - 1) // 6 or not ensembles * walkers <= unique <= samples
                or not source_size or extra_cosmology not in (0, 6)):
            raise ValueError(f"Invalid MCMC chain archive header: {path}")
        source = os.fsdecode(read_exact(stream, source_size))
        title = read_exact(stream, title_size).decode("utf-8")
        mask = map_array(stream, "u1", (ensembles, walkers, (steps + 7) // 8))
        pad_numeric_block(stream)
        cosmology = map_array(stream, "<f8", (extra_cosmology,))
        diagnostics = map_array(stream, "<f8", (ensembles, 6))
        states = map_array(stream, "<f8", (6, unique))
        expect_end(stream)
    return ChainArchive(fit_id, source_id, source, title, (ensembles, steps, walkers),
                        burn, seed, mask, diagnostics, states, cosmology)


def read_curves(path, chain):
    with path.open("rb") as stream:
        magic, version, fit_id, radius_count, added, draws, unique, best = CURVE_HEADER.unpack(
            read_exact(stream, CURVE_HEADER.size)
        )
        if (magic != CURVE_MAGIC or version != VERSION or not 2 <= radius_count <= 436
                or added > radius_count or not 1 <= draws <= math.prod(chain.shape)
                or not 1 <= unique <= draws + 1 or best >= unique):
            raise ValueError(f"Invalid ESD curve archive header: {path}")
        if fit_id != chain.fit_id:
            raise ValueError("Chain and ESD curve archives have different fit IDs")
        pad_numeric_block(stream)
        extra_radii = map_array(stream, "<f8", (added,))
        indices = map_array(stream, "<u4", (draws,))
        rows = map_array(stream, "<u4", (draws,))
        components = map_array(stream, "<f8", (2, radius_count, unique))
        expect_end(stream)
    if (indices[-1] >= math.prod(chain.shape) or np.any(indices[1:] <= indices[:-1])
            or np.any(rows >= unique)):
        raise ValueError("Invalid ESD curve sample indices")
    return CurveArchive(radius_count, extra_radii, indices, rows, best, components)


@contextmanager
def lock_directory(path, mode):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        fcntl.flock(fd, mode)
        yield
    finally:
        os.close(fd)


def archive_names(stem):
    return f"{stem}_chain.bin", f"{stem}_profile_draws.bin"


def load_fit(archive):
    archive = Path(archive).absolute()
    with lock_directory(archive.parent, fcntl.LOCK_SH):
        chain_name, curve_name = archive_names(archive.name)
        chain = read_chain(archive / chain_name)
        curves = read_curves(archive / curve_name, chain)
    if chain.fit_id != archive_id(chain, curves):
        raise ValueError("Fit archive checksum mismatch")
    chain.measurement_path = str(archive.parent / chain.measurement_path)
    return chain, curves


def write_chain(path, chain):
    source = os.fsencode(chain.measurement_path)
    title = chain.title.encode("utf-8")
    with path.open("wb") as stream:
        stream.write(CHAIN_HEADER.pack(CHAIN_MAGIC, VERSION, chain.fit_id, chain.measurement_id,
                                      *chain.shape, chain.burn_in, chain.states.shape[-1], chain.seed,
                                      len(source), len(title), len(chain.legacy_cosmology)))
        stream.write(source)
        stream.write(title)
        chain.mask.tofile(stream)
        pad_numeric_block(stream, writing=True)
        chain.legacy_cosmology.tofile(stream)
        chain.diagnostics.tofile(stream)
        chain.states.tofile(stream)


def write_curves(path, curves, fit_id):
    with path.open("wb") as stream:
        stream.write(CURVE_HEADER.pack(CURVE_MAGIC, VERSION, fit_id, curves.radius_count,
                                       len(curves.extra_radii), len(curves.chain_indices),
                                       curves.components.shape[-1], curves.best_row))
        pad_numeric_block(stream, writing=True)
        curves.extra_radii.tofile(stream)
        curves.chain_indices.tofile(stream)
        curves.rows.tofile(stream)
        curves.components.tofile(stream)


def remove_generation(directory, names):
    for name in (*names, ".publish"):
        (directory / name).unlink(missing_ok=True)
    directory.rmdir()


def save_fit(archive, samples, diagnostics, curves, context):
    archive = Path(archive).absolute()
    archive.parent.mkdir(parents=True, exist_ok=True)
    mask, states = encode_chain(samples)
    diagnostics = np.asarray(diagnostics, dtype="<f8", order="C")
    cosmology = np.asarray(context.get("legacy_cosmology", ()), dtype="<f8")
    shape = samples.shape[:-1]
    source_id, burn, seed = (context[key] for key in ("measurement_id", "burn_in", "seed"))
    source = os.path.relpath(context["measurement_path"], archive.parent)
    chain = ChainArchive(b"", source_id, source, context["title"], shape, burn, seed, mask,
                         diagnostics, states, cosmology)
    chain.fit_id = archive_id(chain, curves)

    generation_prefix = f".{archive.name}.fit-"
    generation = Path(mkdtemp(prefix=generation_prefix, dir=archive.parent))
    published = False
    names = archive_names(archive.name)
    try:
        # Close both files before changing the current fit link.
        write_chain(generation / names[0], chain)
        write_curves(generation / names[1], curves, chain.fit_id)
        link = generation / ".publish"
        link.symlink_to(generation.name, target_is_directory=True)
        with lock_directory(archive.parent, fcntl.LOCK_EX):
            previous = None
            if archive.is_symlink():
                target = Path(os.readlink(archive))
                if target.name == str(target) and target.name.startswith(generation_prefix):
                    previous = archive.parent / target
            elif archive.exists():
                raise FileExistsError(f"Fit archive path already exists: {archive}")
            os.replace(link, archive)
            published = True
            if previous is not None and previous.is_dir() and not previous.is_symlink():
                try:
                    remove_generation(previous, names)
                except OSError:
                    pass
    finally:
        if not published:
            remove_generation(generation, names)
