"""Minimal PicoQuant PTU reader (header tags + T3 records) for the Luminosa PDA-23."""
import struct
import numpy as np

_TY = {
    0xFFFF0008: "Empty8", 0x00000008: "Bool8", 0x10000008: "Int8",
    0x11000008: "BitSet64", 0x12000008: "Color8", 0x20000008: "Float8",
    0x21000008: "TDateTime", 0x2001FFFF: "Float8Array", 0x4001FFFF: "AnsiString",
    0x4002FFFF: "WideString", 0xFFFFFFFF: "BinaryBlob",
}


def read_header(path):
    """Return (tags dict, byte offset of the first record)."""
    tags = {}
    with open(path, "rb") as f:
        magic = f.read(8).rstrip(b"\0").decode()
        if magic != "PQTTTR":
            raise ValueError(f"{path}: not a PTU file ({magic!r})")
        f.read(8)  # version
        while True:
            ident = f.read(32).rstrip(b"\0").decode(errors="replace")
            idx, typ = struct.unpack("<iI", f.read(8))
            raw = f.read(8)
            name = ident if idx < 0 else f"{ident}({idx})"
            kind = _TY.get(typ, hex(typ))
            if kind in ("AnsiString", "WideString", "BinaryBlob", "Float8Array"):
                n = struct.unpack("<q", raw)[0]
                blob = f.read(n)
                if kind == "AnsiString":
                    val = blob.rstrip(b"\0").decode(errors="replace")
                elif kind == "WideString":
                    val = blob.decode("utf-16-le", errors="replace").rstrip("\0")
                elif kind == "Float8Array":
                    val = np.frombuffer(blob, "<f8")
                else:
                    val = blob
            elif kind in ("Float8", "TDateTime"):
                val = struct.unpack("<d", raw)[0]
            elif kind == "Bool8":
                val = bool(struct.unpack("<q", raw)[0])
            else:
                val = struct.unpack("<q", raw)[0]
            tags[name] = val
            if ident == "Header_End":
                return tags, f.tell()


def read_t3(path, max_records=None):
    """Decode generic (HydraHarp2/MultiHarp/PicoHarp330) T3 records.

    Returns dict with absolute sync number `nsync` (int64), `chan` (uint8,
    routing channel as stored, 0-based), `dtime` (uint16 TCSPC bin) and
    the marker records separately.
    """
    tags, off = read_header(path)
    n = tags["TTResult_NumberOfRecords"]
    if max_records is not None:
        n = min(n, max_records)
    rec = np.fromfile(path, dtype="<u4", count=n, offset=off)
    special = (rec >> 31).astype(np.uint8)
    chan = ((rec >> 25) & 0x3F).astype(np.uint8)
    dtime = ((rec >> 10) & 0x7FFF).astype(np.uint16)
    nsync = (rec & 0x3FF).astype(np.int64)
    ovfl = (special == 1) & (chan == 0x3F)
    # overflow records carry the number of overflows in nsync (0 means 1 in old format)
    ovcount = np.where(ovfl, np.where(nsync == 0, 1, nsync), 0)
    offset = np.cumsum(ovcount) * 1024
    nsync_abs = nsync + offset
    photon = special == 0
    marker = (special == 1) & (chan != 0x3F) & (chan > 0)
    return {
        "tags": tags,
        "nsync": nsync_abs[photon],
        "chan": chan[photon],
        "dtime": dtime[photon],
        "marker_nsync": nsync_abs[marker],
        "marker_chan": chan[marker],
    }
