import os
import subprocess
import numpy as np


def sloane(n=None):
    """Return RA and DEC of point coverings of sphere (Hardin, Sloane & Smith 1994)

    Parameters:
    ----------

    n : int, np.int32
        number of points in covering

    Returns:
    -------

    ra : np.float64
        right ascension (deg)
    dec : np.float64
        declination (deg)

    Comments:
    --------

    Only a limited set of coverings are available. If "n" is outside this
    set, returns None.
    """

    icover_path = os.path.join(os.getenv('OBSERVESIM_DIR'),
                               'bin', 'icover')
    out = subprocess.run([icover_path, str(n)], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    if(out.stderr != b''):
        return(None)
    xyz = np.array(out.stdout.splitlines(), dtype=np.float64).reshape(n, 3)
    c = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2)
    ra_rad = (np.arctan2(xyz[:, 1], xyz[:, 0]) + np.pi * 2.) % (np.pi * 2)
    dec_rad = (0.5 * np.pi - np.arctan2(c, xyz[:, 2]))
    ra = ra_rad * 180. / np.pi
    dec = dec_rad * 180. / np.pi

    return(ra, dec)
