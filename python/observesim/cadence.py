import numpy as np
from ortools.constraint_solver import pywrapcp


class Cadence(object):
    """Description of a cadence

    Parameters:
    ----------

    nexposures : np.int32
        number of exposures (default 1)

    lunation : ndarray of np.float32
        maximum lunation for each exposure (default [1.])

    epoch : ndarray of np.float32
        day for exposure (default [0.])

    softness : ndarray of np.float32
        allowance for variation from cadence (default [1.])

    Attributes:
    ----------

    nexposures : np.int32
        number of exposures (default 1)

    lunation : ndarray of np.float32
        maximum lunation for each exposure (default [1.])

    epoch : ndarray of np.float32
        day for exposure (default [0.])

    softness : ndarray of np.float32
        allowance for variation from cadence (default [1.])

"""
    def __init__(self, nexposures=1, lunation=1., epoch=0., softness=1.):
        self.nexposures = np.int32(nexposures)
        self.lunation = np.zeros(self.nexposures, dtype=np.float32) + lunation
        self.epoch = np.zeros(self.nexposures, dtype=np.float32) + epoch
        self.softness = np.zeros(self.nexposures, dtype=np.float32) + softness
        return


class CadenceList(object):
    """List of cadences available

    Parameters:
    ----------

    Attributes:
    ----------

    ncadences : np.int32, int
         number of different cadences

    cadences : list
         list of Cadence objects
"""
    def __init__(self):
        self.ncadences = 0
        self.cadences = []
        return

    def add_cadence(self, *args, **kwargs):
        """Add a cadence to the list of cadences

        Parameters:
        ----------

        nexposures : np.int32
            number of exposures (default 1)

        lunation : ndarray of np.float32
            maximum lunation for each exposure (default [1.])

        epoch : ndarray of np.float32
            day for exposure (default [0.])

        softness : ndarray of np.float32
            allowance for variation from cadence (default [1.])
"""
        cadence = Cadence(*args, **kwargs)
        self.cadences.append(cadence)

    def exposure_consistency(self, epoch_i=None, lunation_i=None,
                             softness_i=None, epoch_j=None,
                             lunation_j=None, softness_j=None):
        """Is exposure i satisfied within constraints of exposure j?

        Parameters:
        ----------

        epoch_i : np.float32
             epoch (in days) of exposure i

        lunation_i : np.float32
             lunation requirement of exposure i

        softness_i : np.float32
             softness of epoch requirement of exposure i

        epoch_j : np.float32
             epoch (in days) of exposure j

        lunation_j : np.float32
             lunation requirement of exposure j

        softness_j : np.float32
             softness of epoch requirement of exposure j
"""
        if(lunation_i < lunation_j):
            return(0)
        if(softness_i > 9999.):
            return(1)
        if(epoch_j - softness_j < epoch_i - softness_i):
            return(0)
        if(epoch_j + softness_j > epoch_i + softness_i):
            return(0)
        return(1)

    def solve_ecm(self, return_solutions=False, ecm=np.array([[1]])):
        """Solve exposure consistency matrix

        Parameters:
        ----------

        ecm : ndarray of np.int32
            exposure consistency matrix

        return_solutions: boolean
            return list of solutions? (default False)

        Returns:
        -------

        ok : int
            1 if there is a solution, 0 otherwise

        solutions : list (if return_solutions is True)
            list of lists of indices into cadence #2

        Notes:
        -----

        Exposure consistency matrix is matrix between two cadences,
        where E_ij is 1 if exposure i from cadence #1 is consistent
        with exposure j of cadence #2, 0 otherwise.

        A solution is a set of N_1 exposures within cadence 2, which
        satisfy the requirements for cadence 1.

        Finds solutions using constraint programming solver from ortools.
"""
        solver = pywrapcp.Solver("ECM")
        ni, nj = ecm.shape
        if(nj < ni):
            if(return_solutions):
                return(0, [])
            else:
                return(0)
        icadence = list()
        for i in np.arange(ni):
            j = [int(x) for x in np.where(ecm[i, :] > 0)[0]]
            if(len(j) == 0):
                if(return_solutions):
                    return(0, [])
                else:
                    return(0)
            icadence.append(solver.IntVar(j, "i={i}".format(i=str(i))))
        solver.Add(solver.AllDifferent(icadence))
        db = solver.Phase(icadence, solver.CHOOSE_FIRST_UNBOUND,
                          solver.ASSIGN_MIN_VALUE)
        solver.Solve(db)
        count = 0
        if(return_solutions):
            solutions = list()
        while solver.NextSolution():
            count += 1
            if(return_solutions):
                solution = list()
                for ii in icadence:
                    solution.append(ii.Value())
                solutions.append(solution)
        if(return_solutions):
            return(count, solutions)
        else:
            return(count)

    def ecm(self, one, two, jstart=0):
        """Create exposure consistency matrix of cadence #1 within cadence #2

        Parameters:
        ----------

        one : int
            index for cadence #1

        two : int
            index for cadence #2

        jstart : int
            which exposure in cadence #2 to start at

        Notes:
        -----

        Exposure consistency matrix is matrix between two cadences,
        where E_ij is 1 if exposure i from cadence #1 is consistent
        with exposure j of cadence #2, 0 otherwise.

        jstart chooses which exposure of cadence #2 to align the first
        exposure of cadence #1 with.
"""
        ecm = np.zeros((self.cadences[one].nexposures,
                        self.cadences[two].nexposures),
                       dtype=np.int32)
        for ic in np.arange(self.cadences[one].nexposures):
            for jc in np.arange(self.cadences[two].nexposures):
                ecm[ic, jc] = self.exposure_consistency(epoch_i=self.cadences[one].epoch[ic],
                                                        lunation_i=self.cadences[one].lunation[ic],
                                                        softness_i=self.cadences[one].softness[ic],
                                                        epoch_j=self.cadences[two].epoch[jc] - self.cadences[two].epoch[jstart],
                                                        lunation_j=self.cadences[two].lunation[jc],
                                                        softness_j=self.cadences[two].softness[jc])
        return(ecm)

    def cadence_consistency(self, one, two, return_solutions=False):
        """Is cadence #1 consistent with cadence #2?

        Parameters:
        ----------

        one : int
            index for cadence #1

        two : int
            index for cadence #2

        return_solutions: boolean
            return list of solutions? (default False)

        Returns:
        -------

        ok : int
            1 if there is a solution, 0 otherwise

        solutions : list (if return_solutions is True)
            list of lists of indices into cadence #2

        Notes:
        -----

        A solution is a set of N_1 exposures within cadence 2, which
        satisfy the requirements for cadence 1.

"""
        count = 0
        solutions = []

        for jstart in np.arange(self.cadences[two].nexposures):
            ecm = self.ecm(one, two, jstart=jstart)
            if(return_solutions):
                (tmp_count, tmp_solutions) = self.solve_ecm(ecm=ecm,
                                                            return_solutions=return_solutions)
                tmp_solutions = tmp_solutions
                if(tmp_count > 0):
                    count = count + tmp_count
                    solutions = solutions + tmp_solutions
            else:
                tmp_count = self.solve_ecm(ecm=ecm,
                                           return_solutions=return_solutions)
                count = count + tmp_count

        if(return_solutions):
            return(count, solutions)
        else:
            return(count)
