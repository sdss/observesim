import os

import yaml
import numpy as np
import fitsio

"""Scheduler module class.

Dependencies:

 numpy

"""


class SN2(object):
    """Class to incapsulate SN2 fiducials and report appropriately random SN2.

       A cumulative distribution function (CDF) constructed from historical data
       and a grid of the same length with values from [0, 1] are required.

       When called, SN2 will report a value appropriately sampled from
       the historical distribution.
    """

    def __init__(self, fid_cdf=None, fid_grid=None):
        """
         Parameters:
        ----------
            fid_cdf : np.float64, array
                The fiducial cdf created from hisotorical data

            fid_grid : np.float64, array
                A grid of len(fid_cdf) with values from [0,1] to allow proper sampling of the CDF
                using a random variable.

        """
        assert fid_cdf is not None, "must init with fiducial CDF!"
        assert fid_grid is not None, "must init with fiducial grid!"

        assert len(fid_cdf) == len(fid_grid), "fiducial CDF and grid must be the same length!"

        self._fid_cdf = fid_cdf
        self._fid_grid = fid_grid

    def _sampleCDF(self, mu):
        """Pull a sample from self.fid_cdf for mu [0,1]

        Parameters:
        ----------

        mu : float
            A random number drawn from [0,1].

        Returns:
        -------

        SN2 : np.float64
            SN2 correpsonding to the input random variable

        """
        idx = np.argmin(np.abs(mu-self._fid_grid))
        return self._fid_cdf[idx]

    def __call__(self, mu=None):
        """obtain a SN2 from fid_cdf

        Parameters:
        ----------

        mu : float, optinal
            A random number drawn from [0,1], if not passed a random number will be chosen.
            This flexibility is intentional, so that for example a value closer to 0 can be
            chosen if the weather is bad, or two close values can be chosen for back-to-back
            observations. Such choices should be made elsewhere however.

        Returns:
        -------

        SN2 : np.float64, array
            An array of SN2 correpsonding to the input random or psuedo-random variable(s)

        """

        if mu is None:
            mu = np.random.rand()
        elif np.isscalar(mu):
            mu = np.array([mu])

        return np.array([self._sampleCDF(x) for x in mu])


class simple_SN2(object):
    """Class to predict SN2 based on airmass
    """

    def __init__(self,
                 b_fit_coef=None, b_noise_coef=None, b_cloudy_shift=0,
                 r_fit_coef=None, r_noise_coef=None, r_cloudy_shift=0,
                 ap_fit_coef=None, ap_noise_coef=None, ap_cloudy_shift=0,
                 ):
        """
         Parameters:
        ----------
            fit_coef : list
                Coefficients for N-degree 1-D polynomial fit for S/N vs airmass

            noise_coef : list
                Coefficients for N-degree 1-D polynomial fit to noise model

        """

        self.b_model = np.poly1d(b_fit_coef)
        self.b_noise = np.poly1d(b_noise_coef)
        self.b_cloudy_shift = b_cloudy_shift
        self.r_model = np.poly1d(r_fit_coef)
        self.r_noise = np.poly1d(r_noise_coef)
        self.r_cloudy_shift = r_cloudy_shift
        self.ap_model = np.poly1d(ap_fit_coef)
        self.ap_noise = np.poly1d(ap_noise_coef)
        self.ap_cloudy_shift = ap_cloudy_shift

    def __call__(self, am=1.0, cloudy=False):
        """obtain a SN2 from fid_cdf

        Parameters:
        ----------

        am : float
            airmass
        cloudy: bool
            apply a negative shift for cloudy weather

        Returns:
        -------

        SN2 : np.float64, array
            An array of SN2 correpsonding to the input airmass

        """

        if np.isscalar(am):
            am = np.array([am])
        if cloudy:
            b_shift = -1.5*self.b_cloudy_shift
            r_shift = -1.5*self.r_cloudy_shift
            ap_shift = -1.5*self.ap_cloudy_shift
        else:
            b_shift = self.b_cloudy_shift
            r_shift = self.r_cloudy_shift
            ap_shift = self.ap_cloudy_shift
        N = len(am)

        random_shared = np.random.randn(N)
        # random noise could push it negative, can't have that!
        b_raw = self.b_model(am) + random_shared*self.b_noise(am) + b_shift
        r_raw = self.r_model(am) + random_shared*self.r_noise(am) + r_shift
        ap_raw = self.ap_model(am) + random_shared*self.ap_noise(am) + ap_shift

        return(np.clip(r_raw, 0, None), 
               np.clip(b_raw, 0, None), 
               np.clip(ap_raw, 0, None))


class Observe(object):
    """Observe class.

    Used to simulate observations of tiles.

    Methods:
    -------

    result() : return result of an observation

    """
    def __init__(self, defaultExp=None, loc="apo"):
        """Return simulated result of an observation

        Parameters:
        ----------

        defaultExp : float
            default exposure length in days

        cadencelist :
            the cadencelist object from scheduler

        cadences :
            the cadences from scheduler (of length fields, matching that index)
        """

        self._result_dtype = [('apgSN2', np.float32),
                              ('rSN2', np.float32),
                              ('bSN2', np.float32),
                              ('mjd', np.float64),
                              ('duration', np.float64)]
        self._result0 = np.zeros(1, dtype=self._result_dtype)

        if defaultExp is None:
            self.defaultExp = np.float32(15. / 60. / 24.)
        else:
            self.defaultExp = defaultExp
        # self.defaultSN2 = 3000.

        # fid_file = '/'.join(os.path.realpath(__file__).split('/')[0:-3]) + "/data/sn_fiducials.fits"

        # fid_data = fitsio.read(fid_file)

        # self.apgSN2 = SN2(fid_cdf=fid_data["apg_fid"], fid_grid=fid_data["apg_grid"])
        # self.rSN2 = SN2(fid_cdf=fid_data["r_fid"], fid_grid=fid_data["r_grid"])
        # self.bSN2 = SN2(fid_cdf=fid_data["b_fid"], fid_grid=fid_data["b_grid"])

        model_file = '/'.join(os.path.realpath(__file__).split('/')[0:-1]) + f"/etc/sn_models_{loc.lower()}.yml"
        models = yaml.load(open(model_file))

        # self.SN2 = simple_SN2(b_fit_coef=models["b"],
        #                       b_noise_coef=models["bnoise"],
        #                       b_cloudy_shift=models["bshift"],
        #                       r_fit_coef=models["r"],
        #                       r_noise_coef=models["rnoise"],
        #                       r_cloudy_shift=models["rshift"],
        #                       ap_fit_coef=models["ap"],
        #                       ap_noise_coef=models["apnoise"],
        #                       ap_cloudy_shift=models["apshift"])
        self.SN2 = simple_SN2(**models)

    def _amSN(self, am, cloudy=False):
        r, b, apg = self.SN2(am, cloudy=cloudy)
        return r, b, 2*apg  # 2x apg exp for each boss exp

    def result(self, field_pk=None, mjd=None, airmass=1,
               epochidx=None, cloudy=False, **kwargs):
        """Return simulated result of an observation

        Parameters:
        ----------

        field_pk : int, np.int32
            pk of field to observe

        mjd : float, np.float64
            MJD start time of observation (days)

        airmass: float, np.float64
            The airmass of the observation

        **kwargs
            catch old duration kwarg quietly?

        Returns:
        -------

        result : dict
            has keys 'sn2', 'mjd', 'duration'
        """

        r, b, apg = self._amSN(airmass, cloudy=cloudy)

        fresult = self._result0
        fresult['apgSN2'] = apg
        fresult['rSN2'] = r
        fresult['bSN2'] = b
        fresult['mjd'] = mjd
        fresult['duration'] = self.defaultExp
        return(fresult)
