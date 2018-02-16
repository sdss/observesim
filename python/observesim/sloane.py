import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import PyAstronomy.pyasl as pyasl
import mpl_toolkits.basemap as basemap


def _arrayify(quantity=None):
    """Cast quantity as ndarray of numpy.float64"""
    try:
        length = len(quantity)
    except TypeError:
        length = 1
    return np.zeros(length, dtype=np.float64) + quantity


def _lb2radec(l=None, b=None):
    """Convert Galactic to equatorial"""
    l = _arrayify(l)
    b = _arrayify(b)
    ra = np.zeros(len(l), dtype=np.float64)
    dec = np.zeros(len(l), dtype=np.float64)
    for i in np.arange(len(l)):
        tany_num = (np.sin(np.deg2rad(l[i] - 123.)))
        tany_den = (np.cos(np.deg2rad(l[i] - 123.)) *
                    np.sin(np.deg2rad(27.4)) -
                    np.tan(np.deg2rad(b[i])) *
                    np.cos(np.deg2rad(27.4)))
        y = np.arctan2(tany_num, tany_den)
        ra1950 = np.rad2deg(y) + 12.25
        sind = (np.sin(np.deg2rad(b[i])) * np.sin(np.deg2rad(27.4)) +
                np.cos(np.deg2rad(b[i])) * np.cos(np.deg2rad(27.4)) *
                np.cos(np.deg2rad(l[i] - 123.)))
        dec1950 = np.rad2deg(np.arcsin(sind))
        ra[i], dec[i] = pyasl.astroTimeLegacy.precess(ra1950, dec1950,
                                                      1950., 2000.)
        ra[i] = ra[i] % 360.
    return (ra, dec)


def _radec2lb(ra=None, dec=None):
    """Convert equatorial to Galactic"""
    ra = _arrayify(ra)
    dec = _arrayify(dec)
    l = np.zeros(len(ra), dtype=np.float64)
    b = np.zeros(len(ra), dtype=np.float64)
    for i in np.arange(len(ra)):
        ra1950, dec1950 = pyasl.astroTimeLegacy.precess(ra[i], dec[i],
                                                        2000., 1950.)
        tanx_num = (np.sin(np.deg2rad(192.25 - ra1950)))
        tanx_den = (np.cos(np.deg2rad(192.25 - ra1950)) *
                    np.sin(np.deg2rad(27.4)) -
                    np.tan(np.deg2rad(dec1950)) *
                    np.cos(np.deg2rad(27.4)))
        l[i] = 303. - np.rad2deg(np.arctan2(tanx_num, tanx_den))
        l[i] = l[i] % 360.
        sinb = (np.sin(np.deg2rad(dec1950)) * np.sin(np.deg2rad(27.4)) +
                np.cos(np.deg2rad(dec1950)) * np.cos(np.deg2rad(27.4)) *
                np.cos(np.deg2rad(192.25 - ra1950)))
        b[i] = np.rad2deg(np.arcsin(sinb))
    return (l, b)


class Sloane(object):
    """RA and DEC of point coverings of sphere (Hardin, Sloane & Smith 1994)

    Parameters:
    ----------

    n : int, np.int32
        number of points in covering

    Attributes:
    ----------

    ra : np.float64
        right ascension (deg)
    dec : np.float64
        declination (deg)

    Methods:
    -------

    plot() : plot ra and dec
    squash() : squash in Galactic latitude

    Comments:
    --------

    Only a limited set of coverings are available. If "n" is outside this
    set, ra and dec are set to None.
"""
    def __init__(self, n=None, alignment='Galactic', radius=1.5):
        self.n = n
        self.alignment = alignment
        (lon, lat) = self._uniform(n=self.n)
        if(lon is None):
            self.ra = None
            self.dec = None
            return
        if(self.alignment == 'Galactic'):
            (ra, dec) = _lb2radec(l=lon, b=lat)
        if(self.alignment == 'Equatorial'):
            (ra, dec) = (lon, lat)
        self.ra = ra
        self.dec = dec
        self.radius = radius
        return

    def _uniform(self, n=None):
        icover_path = os.path.join(os.getenv('OBSERVESIM_DIR'),
                                   'bin', 'icover')
        out = subprocess.run([icover_path, str(n)], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        if(out.stderr != b''):
            return(None, None)
        xyz = np.array(out.stdout.splitlines(), dtype=np.float64).reshape(n, 3)
        c = np.sqrt(xyz[:, 0]**2 + xyz[:, 1]**2)
        ra_rad = (np.arctan2(xyz[:, 1], xyz[:, 0]) + np.pi * 2.) % (np.pi * 2)
        dec_rad = (0.5 * np.pi - np.arctan2(c, xyz[:, 2]))
        ra = ra_rad * 180. / np.pi
        dec = dec_rad * 180. / np.pi
        print(len(ra))

        return(ra, dec)

    def _squash_func(self, l=None, b=None):
        l = _arrayify(l)
        b = _arrayify(b)
        bsign = np.ones(len(b))
        ineg = np.where(b < 0.)[0]
        bsign[ineg] = -1.
        babs = np.abs(b)
        newbabs = ((babs / 90.)**2) * 90.
        newb = bsign * newbabs
        newl = l
        return(newl, newb)

    def squash(self, func=None):
        """Squash tiling in Galactic coordinates

        Parameters:
        ----------

        func : function
            function taking (l, b), returning (newl, newb)

        Comments:
        --------

        By default, squashes just in b to make denser tiling
        in Galactic plane.
"""
        l, b = _radec2lb(ra=self.ra, dec=self.dec)
        if(func is None):
            (newl, newb) = self._squash_func(l=l, b=b)
        else:
            (newl, newb) = func(l=l, b=b)
        ra, dec = _lb2radec(l=newl, b=newb)
        self.ra = ra
        self.dec = dec
        return

    def deccut(self, dec=None, above=True):
        if(above is True):
            ikeep = np.where(self.dec > dec)[0]
        else:
            ikeep = np.where(self.dec < dec)[0]
        self.ra = self.ra[ikeep]
        self.dec = self.dec[ikeep]

    def _convert_radec(self, m, ra, dec):
        return m(((360. - ra) + 180.) % 360., dec, inverse=False)

    def plot(self):
        """Plot ra and dec"""
        m = basemap.Basemap(projection='moll', lon_0=270, resolution='c')

        # draw parallels and meridians.
        m.drawparallels(np.arange(-90., 120., 30.),
                        linewidth=0.5,
                        labels=[1, 0, 0, 0],
                        labelstyle='+/-')
        m.drawmeridians(np.arange(0., 420., 60.), linewidth=0.5)
        m.drawmapboundary()

        boundary = 90.

        ncirc = 60

        theta = np.arange(ncirc) * np.pi * 2. / np.float32(ncirc - 1)

        for (ra, dec) in zip(self.ra, self.dec):
            cra = (ra + self.radius / np.cos(dec * np.pi / 180.) *
                   np.cos(theta))
            cdec = dec + self.radius * np.sin(theta)

            nbelow = np.count_nonzero(cra < boundary)
            if(nbelow > 0):
                ibelow = np.nonzero(cra < boundary)[0]
                (xx, yy) = self._convert_radec(m, cra[ibelow], cdec[ibelow])
                plt.plot(xx, yy, linewidth=0.25, color='blue')

            nabove = np.count_nonzero(cra > boundary)
            if(nabove > 0):
                iabove = np.nonzero(cra > boundary)[0]
                (xx, yy) = self._convert_radec(m, cra[iabove], cdec[iabove])
                plt.plot(xx, yy, linewidth=0.25, color='blue')
