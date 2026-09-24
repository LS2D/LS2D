#
# This file is part of LS2D.
#
# Copyright (c) 2017-2026 Wageningen University & Research
# Author: Bart van Stratum (WUR)
#
# LS2D is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LS2D is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with LS2D.  If not, see <http://www.gnu.org/licenses/>.
#

"""
Partial reads of small lat/lon boxes from the Google ARCO-ERA5 Zarr stores.

The ARCO stores are chunked as global fields, e.g. (time=1, level=18, lat=721, lon=1440).
Reading a small box through xarray/zarr downloads and decompresses all of those.
However, every chunk is a Blosc frame of independently compressed ~512 kB blocks,
with a table of block offsets in its header. So instead we:
  1. range-read the Blosc header + block offset table of every chunk we need,
  2. range-read only the blocks that overlap the latitude rows of our box,
  3. decode each block separately by wrapping it in a single-block Blosc frame.
"""

# Python modules
import struct

# Third party modules
import numpy as np
from numcodecs import blosc

# LS2D modules
from ls2d.core.messages import *

# Blosc header is 16 bytes, followed by one int32 block offset per block.
_BLOSC_HEADER = 16
# Range read for header + offset table; 2 kB covers up to 508 blocks.
_HEAD_BYTES = 2048


def get_layout(ds, name):
    """
    Get and check the chunk layout of variable `name` in the (metadata only) ARCO dataset `ds`.
    Returns dict with the number of levels and levels per chunk (both 1 for 2D fields).
    """

    da = ds[name]
    enc = da.encoding
    compressor = enc['compressors'][0] if 'compressors' in enc else enc['compressor']
    chunks = enc['chunks']

    ok = (
        da.dims[0] == 'time'
        and da.dims[-2:] == ('latitude', 'longitude')
        and da.ndim in (3, 4)
        and chunks[0] == 1
        and tuple(chunks[-2:]) == (ds.sizes['latitude'], ds.sizes['longitude'])
        and enc['dtype'] == np.float32
        and compressor.codec_id == 'blosc'
    )
    if not ok:
        error(f'ARCO variable "{name}" has an unsupported layout: dims={da.dims}, chunks={chunks}')

    if da.ndim == 4:
        return dict(nlev=da.shape[1], nlev_chunk=chunks[1])
    else:
        return dict(nlev=1, nlev_chunk=1)


def _decode_block(header, raw, iblock, nblocks):
    """
    Decode one Blosc block by wrapping it into a single-block Blosc frame.
    """
    ver, verlz, flags, typesize, nbytes, blocksize, _ = struct.unpack('<BBBBiii', header)
    leftover = nbytes % blocksize
    bsize = leftover if (iblock == nblocks - 1 and leftover) else blocksize
    frame = struct.pack('<BBBBiiii', ver, verlz, flags, typesize, bsize, blocksize, 20 + len(raw), 20)
    return blosc.decompress(frame + raw)


def read_rows(fs, fields, it, ilat0, nlat, nlat_g, nlon_g, batch_size=512):
    """
    Read latitude rows `ilat0:ilat0+nlat` (all longitudes, all levels) at time index `it`.

    Arguments:
        fs : gcsfs.GCSFileSystem
        fields : list of (store_path, variable name, layout from `get_layout()`)
        it : time index in the Zarr store
        ilat0, nlat : first latitude index and number of latitude rows
        nlat_g, nlon_g : global number of latitudes and longitudes
        batch_size : number of concurrent HTTP range requests
    Returns:
        out : dict {name: np.ndarray(nlev, nlat, nlon_g)}
    """

    row_bytes = nlon_g * 4
    out = {}
    chunks = []
    for store, name, lay in fields:
        out[name] = np.full((lay['nlev'], nlat, nlon_g), np.nan, np.float32)
        if lay['nlev'] == 1:
            chunks.append((name, lay, 0, f'{store}/{name}/{it}.0.0'))
        else:
            for k in range(-(-lay['nlev'] // lay['nlev_chunk'])):
                chunks.append((name, lay, k, f'{store}/{name}/{it}.{k}.0.0'))

    paths = [c[3] for c in chunks]

    # 1. Blosc headers + block offset tables.
    heads = fs.cat_ranges(paths, [0] * len(paths), [_HEAD_BYTES] * len(paths), on_error='return', batch_size=batch_size)

    # 2. Figure out which byte ranges we need from every chunk.
    req_paths, req_starts, req_ends, req_meta = [], [], [], []
    for (name, lay, k, path), head in zip(chunks, heads):
        if isinstance(head, FileNotFoundError):
            continue  # Missing chunk == fill value (NaN).
        if isinstance(head, Exception):
            raise head

        header = head[:_BLOSC_HEADER]
        _, _, flags, typesize, nbytes, blocksize, cbytes = struct.unpack('<BBBBiii', header)
        if typesize != 4 or nbytes != lay['nlev_chunk'] * nlat_g * row_bytes:
            error(f'Unexpected Blosc header in {path}: typesize={typesize}, nbytes={nbytes}')

        memcpyed = flags & 0x2
        nblocks = -(-nbytes // blocksize)
        if not memcpyed:
            if _BLOSC_HEADER + 4 * nblocks > len(head):
                error(f'Blosc block table of {path} larger than {_HEAD_BYTES} bytes')
            bstarts = np.frombuffer(head, '<i4', count=nblocks, offset=_BLOSC_HEADER)
            bends = np.append(bstarts[1:], cbytes)

        for lk in range(min(lay['nlev_chunk'], lay['nlev'] - k * lay['nlev_chunk'])):
            s = (lk * nlat_g + ilat0) * row_bytes
            e = s + nlat * row_bytes
            meta = (name, k * lay['nlev_chunk'] + lk, header, s, e, nblocks)
            if memcpyed:
                req_starts.append(_BLOSC_HEADER + s)
                req_ends.append(_BLOSC_HEADER + e)
                req_meta.append(meta + (None,))
            else:
                b0, b1 = s // blocksize, (e - 1) // blocksize
                req_starts.append(int(bstarts[b0]))
                req_ends.append(int(bends[b1]))
                blocks = (b0, b1, bstarts[b0 : b1 + 1] - bstarts[b0], bends[b0 : b1 + 1] - bstarts[b0], blocksize)
                req_meta.append(meta + (blocks,))
            req_paths.append(path)

    # 3. Fetch only the needed blocks, decode, and cut out the latitude rows.
    datas = fs.cat_ranges(req_paths, req_starts, req_ends, on_error='raise', batch_size=batch_size)
    for data, (name, lev, header, s, e, nblocks, blocks) in zip(datas, req_meta):
        if blocks is None:
            buf = data
        else:
            b0, b1, rs, re, blocksize = blocks
            buf = b''.join(_decode_block(header, data[rs[i] : re[i]], b0 + i, nblocks) for i in range(b1 - b0 + 1))
            buf = buf[s - b0 * blocksize : e - b0 * blocksize]
        out[name][lev] = np.frombuffer(buf, np.float32).reshape(nlat, nlon_g)

    return out
