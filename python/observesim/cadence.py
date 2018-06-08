import numpy as np
import itertools
from ortools.constraint_solver import pywrapcp


class Cadence(object):
    """Description of a cadence

    Parameters:
    ----------

    nexposures : np.int32
        number of exposures (default 1)

    lunation : ndarray of np.float32
        maximum lunation for each exposure (default [1.])

    delta : ndarray of np.float32
        day for exposure (default [0.])

    softness : ndarray of np.float32
        allowance for variation from cadence (default [1.])

    Attributes:
    ----------

    nexposures : np.int32
        number of exposures (default 1)

    lunation : ndarray of np.float32
        maximum lunation for each exposure (default [1.])

    delta : ndarray of np.float32
        day for exposure (default [0.])

    epoch : ndarray of np.float32
        day for exposure (default [0.])

    softness : ndarray of np.float32
        allowance for variation from cadence (default [1.])

    Methods:
    -------

    fits(epoch=epoch) : into which parts of cadence an epoch can belong
"""
    def __init__(self, nexposures=1, lunation=1., delta=0., softness=1.):
        self.nexposures = np.int32(nexposures)
        self.lunation = np.zeros(self.nexposures, dtype=np.float32) + lunation
        self.delta = np.zeros(self.nexposures, dtype=np.float32) + delta
        self.epoch = self.delta.cumsum()
        self.softness = np.zeros(self.nexposures, dtype=np.float32) + softness
        return

    def _arrayify(self, quantity=None, dtype=np.float64):
        """Cast quantity as ndarray of numpy.float64"""
        try:
            length = len(quantity)
        except TypeError:
            length = 1
        return np.zeros(length, dtype=dtype) + quantity

    def fits(self, epoch=None, lunation=-1., check_lunation=True):
        fit = (np.abs(epoch - self.epoch) <= self.softness)
        if(check_lunation):
            fit = fit & (lunation <= self.lunation)
        ifit = np.where(fit)[0]
        return(ifit)

    def evaluate(self, epochs=None, lunations=None, check_lunation=True):
        iobs = np.zeros(self.nexposures) - 1
        icad = np.zeros(len(epochs)) - 1
        for indx, epoch, lunation in zip(np.arange(len(epochs)), epochs,
                                         lunations):
            ifit = self.fits(epoch=epoch, lunation=lunation,
                             check_lunation=check_lunation)
            if(len(ifit) > 0):
                iok = np.where(iobs[ifit] == -1)[0]
                if(len(iok) > 0):
                    iobs[ifit[iok[0]]] = indx
                    icad[indx] = ifit[iok[0]]
        return(icad, iobs)

    def evaluate_next(self, mjd_past=None, mjd_next=None,
                      lunation_past=None, lunation_next=None,
                      check_lunation=True):
        mjd_past = self._arrayify(mjd_past)
        icad, iobs = self.evaluate(mjd_past - mjd_past.min(),
                                   lunation_past,
                                   check_lunation=check_lunation)
        ifit = self.fits(epoch=mjd_next - mjd_past.min(),
                         lunation=lunation_next, check_lunation=check_lunation)
        iok = np.where(iobs[ifit] == -1)[0]
        if(len(iok) > 0):
            return True
        else:
            return False


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

        delta : ndarray of np.float32
            day for exposure (default [0.])

        softness : ndarray of np.float32
            allowance for variation from cadence (default [1.])
"""
        cadence = Cadence(*args, **kwargs)
        self.cadences.append(cadence)
        self.ncadences = self.ncadences + 1

    def check_exposures(self, one=None, two=None, indx2=None):
        """Is exposure set in cadence two consistent with cadence one?

        Parameters:
        ----------

        one : int
            index for cadence #1

        two : int
            index for cadence #2

        indx2 : ndarray of np.int32
            exposures in cadence #2
"""
        nexp = len(indx2)

        # Check lunations
        for indx in np.arange(nexp):
            if(self.cadences[one].lunation[indx] <
               self.cadences[two].lunation[indx2[indx]]):
                return(False)

        # Check deltas
        for indx in np.arange(nexp - 1) + 1:
            delta1 = self.cadences[one].delta[indx]
            softness1 = self.cadences[one].softness[indx]
            delta2 = self.cadences[two].delta[indx2[indx - 1] + 1:indx2[indx] + 1].sum()
            softness2 = self.cadences[two].softness[indx2[indx - 1] + 1:indx2[indx] + 1].sum()
            if(delta1 - softness1 > delta2 - softness2):
                return(False)
            if(delta1 + softness1 < delta2 + softness2):
                return(False)

        return(True)

    def cadence_consistency(self, one, two, return_solutions=True):
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
        indx2 = np.arange(self.cadences[two].nexposures)
        sequences = itertools.combinations(indx2,
                                           self.cadences[one].nexposures)
        possibles = []
        for sequence in sequences:
            ok = self.check_exposures(one=one, two=two, indx2=sequence)
            if(ok):
                possibles.append(np.array(sequence))
                if(return_solutions is False):
                    return(True)

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

        target_cadences : ndarray of np.int32
            indices for the target cadences

        field_cadence : int, np.int32
            index for field cadence

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

    def toarray(self):
        """Return cadences as a record array

        Returns:
        -------

        cadences : ndarray
            information on each cadence
"""
        nexps = np.array([c.nexposures for c in self.cadences])
        max_nexp = nexps.max()
        cadence0 = [('cadenceid', np.int32),
                    ('nexposures', np.int32),
                    ('delta', np.float64, max_nexp),
                    ('lunation', np.float32, max_nexp),
                    ('softness', np.float32, max_nexp)]
        cads = np.zeros(self.ncadences, dtype=cadence0)
        for indx, ccadence in zip(np.arange(self.ncadences), self.cadences):
            nexp = self.cadences[indx].nexposures
            cads['cadenceid'][indx] = indx
            cads['nexposures'][indx] = nexp
            cads['delta'][indx, 0:nexp] = self.cadences[indx].delta
            cads['softness'][indx, 0:nexp] = self.cadences[indx].softness
            cads['lunation'][indx, 0:nexp] = self.cadences[indx].lunation
        return(cads)
