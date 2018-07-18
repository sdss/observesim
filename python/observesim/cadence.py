import numpy as np
import itertools
import fitsio
from ortools.constraint_solver import pywrapcp


# Class to define a singleton
class CadenceSingleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(CadenceSingleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Cadence(object):
    """Description of a cadence

    Parameters:
    ----------

    nexposures : np.int32
        number of exposures

    lunation : ndarray of np.float32
        maximum lunation for each exposure

    delta : ndarray of np.float32
        desired offset for each exposure from previous (days)

    delta_min : ndarray of np.float32
        minimum delta to allow (days)

    delta_max : ndarray of np.float32
        maximum delta to allow (days)

    Attributes:
    ----------

    nexposures : np.int32
        number of exposures

    lunation : ndarray of np.float32
        maximum lunation for each exposure

    delta : ndarray of np.float32
        desired offset for each exposure from previous (days)

    delta_min : ndarray of np.float32
        minimum delta to allow (days)

    delta_max : ndarray of np.float32
        maximum delta to allow (days)

    nepochs : np.int32
        number of epochs

    epoch_indx : ndarray of np.int32
        index into delta for first exposure at each epoch

    epoch_delta : ndarray of np.float32
        day for each separate epoch

    epoch_nexposures : ndarray of np.float32
        number of exposures for each separate epoch

    Methods:
    -------
"""
    def __init__(self, nexposures=None, lunation=None, delta=None,
                 delta_min=None, delta_max=None):
        self.nexposures = np.int32(nexposures)
        self.lunation = np.zeros(self.nexposures, dtype=np.float32) + lunation
        self.delta = np.zeros(self.nexposures, dtype=np.float32) + delta
        self.delta_min = (np.zeros(self.nexposures, dtype=np.float32) +
                          delta_min)
        self.delta_max = (np.zeros(self.nexposures, dtype=np.float32) +
                          delta_max)
        self._create_epochs()
        return

    def __str__(self):
        out = "nexposures={nexposures}\n".format(nexposures=self.nexposures)
        out = out + " lunation = "
        for i in np.arange(self.nexposures):
            out = out + " {l:.3f}".format(l=self.lunation[i])
        out = out + "\n"
        out = out + " delta = "
        for i in np.arange(self.nexposures):
            out = out + " {d:.3f}".format(d=self.delta[i])
        out = out + "\n"
        out = out + " delta_min = "
        for i in np.arange(self.nexposures):
            out = out + " {s:.3f}".format(s=self.delta_min[i])
        out = out + "\n"
        out = out + " delta_max = "
        for i in np.arange(self.nexposures):
            out = out + " {s:.3f}".format(s=self.delta_max[i])
        out = out + "\n"
        return(out)

    def _arrayify(self, quantity=None, dtype=np.float64):
        """Cast quantity as ndarray of numpy.float64"""
        try:
            length = len(quantity)
        except TypeError:
            length = 1
        return np.zeros(length, dtype=dtype) + quantity

    def _create_epochs(self):
        epoch_indx = [0]
        self.nepochs = 1
        for indx in np.arange(self.nexposures - 1) + 1:
            if(self.delta[indx] > 0.):
                epoch_indx.append(indx)
                self.nepochs = self.nepochs + 1
        self.epoch_indx = np.array(epoch_indx, dtype=np.int32)
        self.epoch_delta = self.delta[self.epoch_indx]
        self.epoch_delta_min = self.delta_min[self.epoch_indx]
        self.epoch_delta_max = self.delta_max[self.epoch_indx]
        self.epoch_lunation = self.lunation[self.epoch_indx]
        self.epoch_nexposures = np.zeros(self.nepochs, dtype=np.int32)
        for indx in np.arange(self.nepochs - 1):
            self.epoch_nexposures[indx] = (self.epoch_indx[indx + 1] -
                                           self.epoch_indx[indx])
        self.epoch_nexposures[-1] = self.nexposures - self.epoch_indx[-1]
        return

    def evaluate_next(self, mjd_past=None, mjd_next=None,
                      lunation_next=None, check_lunation=True):
        nexposures_past = len(mjd_past)
        if(nexposures_past >= self.nexposures):
            return(False)
        ok_lunation = ((lunation_next < self.lunation[nexposures_past]) |
                       (check_lunation is False))
        if(nexposures_past == 0):
            return(ok_lunation)
        delta = mjd_next - mjd_past[nexposures_past - 1]
        dlo = self.delta_min[nexposures_past]
        dhi = self.delta_max[nexposures_past]
        return(ok_lunation & (delta >= dlo) & (delta <= dhi))


class CadenceList(object, metaclass=CadenceSingleton):
    """List of cadences available (singleton)

    Parameters:
    ----------

    Attributes:
    ----------

    ncadences : np.int32, int
         number of different cadences

    cadences : dictionary
         dictionary of Cadence objects

    Methods:
    -------

    reset() : remove all the current cadences
    add_cadence() : add a new cadence
    check_exposures() : are two exposure sets consistent?
    cadence_consistency(): is cadence #1 consistent with cadence #2?
    pack_targets(): find best way to pack targets into a cadence
    fromarray(): add to cadence list from an ndarray
    fromfits(): add to cadence list from a FITS file
    toarray(): return an ndarray with cadence list

    Notes:
    -----

    This is a singleton, so there can only be one CadenceList defined
    within any session.

"""
    def __init__(self):
        self.ncadences = 0
        self.cadences = dict()
        return

    def reset(self):
        """Reset cadence list to be empty"""
        self.ncadences = 0
        self.cadences = dict()
        return

    def add_cadence(self, name=None, *args, **kwargs):
        """Add a cadence to the list of cadences

        Parameters:
        ----------

        name : string
            dictionary name of cadence

        nexposures : np.int32
            number of exposures (default 1)

        lunation : ndarray of np.float32
            maximum lunation for each exposure (default [1.])

        delta : ndarray of np.float32
            day for exposure (default [0.])

        delta_min : ndarray of np.float32
            allowance for variation from cadence (default [1.])

        delta_max : ndarray of np.float32
            allowance for variation from cadence (default [1.])
"""
        cadence = Cadence(*args, **kwargs)
        self.cadences[name] = cadence
        self.ncadences = self.ncadences + 1

    def check_exposures(self, one=None, two=None, indx2=None, sub1=None):
        """Is exposure set in cadence two consistent with cadence one?

        Parameters:
        ----------

        one : string
            name of cadence #1

        two : string
            name of cadence #2

        indx2 : ndarray of np.int32
            exposures in cadence #2
"""
        eps = 1.e-4  # generously deal with round-off

        if(sub1 is None):
            sub1 = np.arange(self.cadences[one].nexposures)
        nexp = len(indx2)

        # Check lunations
        for indx in np.arange(nexp):
            if(self.cadences[one].lunation[sub1[indx]] <
               self.cadences[two].lunation[indx2[indx]] - eps):
                return(False)

        # Check deltas
        for indx in np.arange(nexp - 1) + 1:
            delta1 = self.cadences[one].delta[sub1[indx]]
            dlo1 = self.cadences[one].delta_min[sub1[indx]]
            dhi1 = self.cadences[one].delta_max[sub1[indx]]
            delta2 = self.cadences[two].delta[indx2[indx - 1] + 1:indx2[indx] + 1].sum()
            dlo2 = self.cadences[two].delta_min[indx2[indx - 1] + 1:indx2[indx] + 1].sum()
            dhi2 = self.cadences[two].delta_max[indx2[indx - 1] + 1:indx2[indx] + 1].sum()
            if(delta1 > 0.):  # normal case
                if(dlo1 >= dlo2 - eps):
                    return(False)
                if(dhi1 <= dhi2 + eps):
                    return(False)
            else:  # adjacent exposures
                if(indx2[indx] > indx2[indx - 1] + 1):  # must be adjacent
                    return(False)
                if(delta2 > 0.):  # must be adjacent
                    return(False)

        return(True)

    def cadence_consistency(self, one, two, return_solutions=True):
        """Is cadence #1 consistent with cadence #2?

        Parameters:
        ----------

        one : string
            name of cadence #1

        two : string
            name of cadence #2

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
        indx2 = np.arange(self.cadences[two].nexposures)

        # Check which exposures you can start on
        possibles = []
        for first in np.arange(self.cadences[two].nexposures -
                               self.cadences[one].nexposures + 1):
            ok = self.check_exposures(one=one, two=two, indx2=[first],
                                      sub1=[0])
            if(ok):
                possibles.append([first])
        if(len(possibles) == 0):
            success = 0
            if(return_solutions):
                return(success, possibles)
            else:
                return(success)

        # Now find sequences starting from there
        for nsub1 in np.arange(self.cadences[one].nexposures - 1) + 2:
            current_possibles = possibles
            possibles = []
            for indx in range(len(current_possibles)):
                possible = current_possibles[indx]
                remaining_start = possible[-1] + 1
                nremaining = self.cadences[two].nexposures - possible[-1] - 1
                ok = 1
                if(nremaining >=
                   self.cadences[one].nexposures - len(possible)):
                    for next_possible in (remaining_start +
                                          np.arange(nremaining)):
                        try_possible = possible.copy()
                        try_possible.append(next_possible)
                        ok = self.check_exposures(one=one, two=two,
                                                  indx2=try_possible,
                                                  sub1=np.arange(nsub1))
                        if(ok):
                            possibles.append(try_possible)
            if(len(possibles) == 0):
                success = 0
                if(return_solutions):
                    return(success, possibles)
                else:
                    return(success)

        success = len(possibles) > 0
        if(return_solutions):
            return(success, possibles)
        else:
            return(success)

    def pack_targets(self, target_cadences=None, field_cadence=None,
                     value=None):
        """Find best way to pack targets into a cadence

        Parameters:
        ----------

        target_cadences : list of strings
            names of the target cadences

        field_cadence : string
            name of field cadence

        value : ndarray of np.float32
            value for each target

        Returns:
        -------

        epoch_targets : ndarray of np.int32
            indices of targets observed for each epoch of field cadence

        Notes:
        -----

        Designed to maximize total "value" of targets observed. Will
        not observe partial cadences.
"""
        ntargets = len(target_cadences)
        nepochs = self.cadences[field_cadence].nexposures
        if(value is None):
            value = np.ones(ntargets)
        else:
            value = np.array(value)

        pack = []
        epochs = [0] * nepochs

        # Find solutions for each target
        for target_cadence in target_cadences:
            count, solns = self.cadence_consistency(target_cadence,
                                                    field_cadence,
                                                    return_solutions=True)
            target = []
            for isoln in solns:
                soln_epochs = epochs.copy()
                for i in isoln:
                    soln_epochs[i] = 1
                    target.append(soln_epochs)
            pack.append(target)

        solver = pywrapcp.Solver("pack_targets")

        # Create variables for each epoch of each solution of each target
        packvar = []
        indxt = 0
        for target in pack:
            indxs = 0
            targetvar = []
            for soln in target:
                indxe = 0
                solnvar = []
                for epoch in soln:
                    name = "{t}-{s}-{e}".format(t=indxt, s=indxt, e=indxe)
                    solnvar.append(solver.BoolVar(name))
                    indxe = indxe + 1
                indxs = indxs + 1
                targetvar.append(solnvar)
            indxt = indxt + 1
            packvar.append(targetvar)

        # Constraint for each solution that all or none of
        # its epochs, and only those epochs, are used
        for target, targetvar in zip(pack, packvar):
            for soln, solnvar in zip(target, targetvar):
                firstvar = None
                for epoch, epochvar in zip(soln, solnvar):
                    if(epoch):
                        if(firstvar):
                            solver.Add(epochvar == firstvar)
                        else:
                            firstvar = epochvar
                    else:
                        solver.Add(epochvar == 0)

        # Constraint for each epoch that only one total
        # target is taken
        for iepoch in range(nepochs):
            e = [solnvar[iepoch] for targetvar in packvar
                 for solnvar in targetvar]
            solver.Add(solver.Sum(e) <= 1)

        # Constraint that only one solution for each target
        # is taken.
        target_valvars = []
        for val, targetvar in zip(value, packvar):
            solnused = []
            for solnvar in targetvar:
                solnused.append(solver.Sum(solnvar) > 0)
            target_got = solver.Sum(solnused)
            target_valvar = target_got * int(val)
            target_valvars.append(target_valvar)
            solver.Add(target_got <= 1)

        allvars = [epochvar for targetvar in packvar
                   for solnvar in targetvar for epochvar in solnvar]

        objective_expr = solver.IntVar(0, int(value.sum()), "value")
        solver.Add(objective_expr == solver.Sum(target_valvars))
        objective = solver.Maximize(objective_expr, 1)

        db = solver.Phase(allvars, solver.CHOOSE_FIRST_UNBOUND,
                          solver.ASSIGN_MIN_VALUE)

        # Create a solution collector.
        collector = solver.LastSolutionCollector()

        # Add the decision variables.
        for allvar in allvars:
            collector.Add(allvar)

        # Add the objective.
        collector.AddObjective(objective_expr)

        success = solver.Solve(db, [objective, collector])
        if(success is False):
            print("Problem in solver.")
            return(None)

        epoch_targets = np.zeros(nepochs, dtype=np.int32) - 1
        if collector.SolutionCount() > 0:
            best_solution = collector.SolutionCount() - 1
            for itarget, targetvar in zip(range(ntargets), packvar):
                for solnvar in targetvar:
                    for iepoch, epochvar in zip(range(nepochs), solnvar):
                        if(collector.Value(best_solution, epochvar)):
                            epoch_targets[iepoch] = itarget

        return(epoch_targets)

    def fromarray(self, cadences_array=None):
        """Add cadences to ccadence list from an array

        Parameters:
        -----------

        cadences_array : ndarray
            ndarray with columns 'nexposures', 'lunation', 'delta',
            'delta_min', 'delta_max', 'cadence'
"""
        for ccadence in cadences_array:
            nexp = ccadence['nexposures']
            self.add_cadence(nexposures=ccadence['nexposures'],
                             lunation=ccadence['lunation'][0:nexp],
                             delta=ccadence['delta'][0:nexp],
                             delta_min=ccadence['delta_min'][0:nexp],
                             delta_max=ccadence['delta_max'][0:nexp],
                             name=ccadence['cadence'].decode().strip())
        return

    def fromfits(self, filename=None):
        """Add cadences to ccadence list from a FITS file

        Parameters:
        -----------

        filename : str
            File name to read from

        Notes:
        -----

        Expects a valid FITS file with columns 'nexposures',
            'lunation', 'delta', 'delta_min', 'delta_max', 'cadence'
"""
        self.cadences_fits = fitsio.read(filename)
        self.fromarray(self.cadences_fits)
        return

    def toarray(self):
        """Return cadences as a record array

        Returns:
        -------

        cadences : ndarray
            information on each cadence
"""
        nexps = np.array([c.nexposures for c in self.cadences.values()])
        max_nexp = nexps.max()
        cadence0 = [('cadence', np.dtype('a20')),
                    ('nexposures', np.int32),
                    ('delta', np.float64, max_nexp),
                    ('lunation', np.float32, max_nexp),
                    ('delta_max', np.float32, max_nexp),
                    ('delta_min', np.float32, max_nexp)]
        cads = np.zeros(self.ncadences, dtype=cadence0)
        names = self.cadences.keys()
        for indx, name in zip(np.arange(self.ncadences), names):
            nexp = self.cadences[name].nexposures
            cads['cadence'][indx] = name
            cads['nexposures'][indx] = nexp
            cads['delta'][indx, 0:nexp] = self.cadences[name].delta
            cads['delta_min'][indx, 0:nexp] = self.cadences[name].delta_min
            cads['delta_max'][indx, 0:nexp] = self.cadences[name].delta_max
            cads['lunation'][indx, 0:nexp] = self.cadences[name].lunation
        return(cads)
