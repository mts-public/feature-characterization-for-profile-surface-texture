import warnings
from pathlib import Path

import numpy as np

# exponents of the supported length units (base: m)
_UNIT_EXPONENTS = {"m": 0, "mm": -3, "um": -6, "nm": -9}


def _unit_factor(unit: str, target: str) -> float:
    """
    conversion factor of a length unit to the target unit
    """
    if unit not in _UNIT_EXPONENTS:
        raise ValueError(f"Unknown unit '{unit}'.")
    return 10.0 ** (_UNIT_EXPONENTS[unit] - _UNIT_EXPONENTS[target])


def read_smd(filepath):
    """
    Reads a profile from a softgauge file (*.smd) according to ISO 5436-2.

    The file consists of four sections separated by ETX (chr(3)):
      1. header with the axis definitions, e.g.
           CX I 8001 mm 1.0E+000 D 5.0E-004  (incremental x-axis, increment)
           CZ A 8001 um 1.0E+000 D           (absolute z-axis)
         fields: axis, axis type, number of points, unit, scale factor,
         data type (and increment for incremental axes)
      2. additional information (e.g. DATE, CREATED-BY)
      3. z-values, one value per line
      4. checksum: sum of all bytes up to the end of the line of the third
         ETX modulo 65535
    Fields can be separated by NUL characters or blanks.

    Parameters
    ----------
        filepath : str or Path
            path of the *.smd file
    Returns
    -------
        z : nd.array, float
            vertical profile values in µm
        L : float
            profile length n*dx in mm
        x : nd.array, float
            x-positions in mm
        dx : float
            step size in x-direction in mm
    """
    data = Path(filepath).read_bytes()

    # sections
    iETX = [i for i, b in enumerate(data) if b == 3]
    if len(iETX) < 3:
        raise ValueError(
            f"'{filepath}' is not a valid smd file (less than 3 ETX separators)."
        )
    header = data[: iETX[0]].replace(b"\x00", b" ").decode("latin-1")
    values = data[iETX[1] + 1 : iETX[2]].decode("latin-1")

    # header: axis definitions
    CX, CZ = None, None
    for line in header.splitlines():
        fields = line.split()
        if fields and fields[0] == "CX":
            CX = fields
        elif fields and fields[0] == "CZ":
            CZ = fields
    if CX is None or CZ is None:
        raise ValueError(
            f"'{filepath}' does not contain a CX and a CZ axis definition."
        )
    if CX[1] != "I" or len(CX) < 7:
        raise ValueError("Only incremental x-axes (CX I ... increment) are supported.")
    if CZ[1] != "A":
        raise ValueError("Only absolute z-axes (CZ A ...) are supported.")
    dx = float(CX[6]) * float(CX[4]) * _unit_factor(CX[3], "mm")
    scale_z = float(CZ[4]) * _unit_factor(CZ[3], "um")

    # z-values
    lines = values.strip().splitlines()
    z = np.empty(len(lines))
    for i, line in enumerate(lines):
        try:
            z[i] = float(line)
        except ValueError:
            raise ValueError(
                f"Invalid z-value in line {i + 1} of the data section."
            ) from None
    z = z * scale_z
    n = z.size
    if n != int(CZ[2]):
        warnings.warn(f"Number of z-values ({n}) differs from the header ({CZ[2]}).")

    # checksum
    if len(iETX) >= 4:
        # end of the line of the third ETX
        iEOL = data.find(b"\n", iETX[2], iETX[3])
        if iEOL < 0:
            iEOL = iETX[2]
        checksum = int(data[iEOL + 1 : iETX[3]].strip())
        if sum(data[: iEOL + 1]) % 65535 != checksum:
            warnings.warn(f"Checksum of '{filepath}' does not match.")

    # x-values
    L = n * dx
    x = np.arange(n) * dx
    return z, L, x, dx
