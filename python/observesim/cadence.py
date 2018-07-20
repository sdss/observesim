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

    instrument : list of str
            instrument for each exposure

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

    instrument : list of str
            instrument for each exposure

    nepochs : np.int32
        number of epochs

    epoch_indx : ndarray of np.int32
        index into delta for first exposure at each epoch

    epoch_nexposures : ndarray of np.float32
        number of exposures for each separate epoch

    Methods:
    -------
"""
    def __init__(self, nexposures=None, lunation=None, delta=None,
                 delta_min=None, delta_max=None, instrument=None):
        self.nexposures = np.int32(nexposures)
        self.lunation = np.zeros(self.nexposures, dtype=np.float32) + lunation
        if(len(instrument) != nexposures):
            print("WHI")
        self.instrument = instrument
        self.delta = np.zeros(self.nexposures, dtype=np.float32) + delta
        self.delta_min = (np.zeros(self.nexposures, dtype=np.float32) +
                          delta_min)
        self.delta_max = (np.zeros(self.nexposures, dtype=np.float32) +
                          delta_max)
        self._create_epochs()
        return

    def __str__(self):
        return(self.exposure_text())

    def exposure_text(self):
        """Display list of exposures as text"""
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
        out = out + " instrument = "
        for i in np.arange(self.nexposures):
            out = out + " {s}".format(s=self.instrument[i])
        out = out + "\n"
        return(out)

    def epoch_text(self):
        """Display list of epochs as text"""
        epoch_indx = self.epoch_indx
        out = "nepochs={nepochs}\n".format(nepochs=self.nepochs)
        out = out + " nexposures = "
        for i in np.arange(self.nepochs):
            out = out + " {l}".format(l=self.epoch_nexposures[i])
        out = out + "\n"
        out = out + " lunation = "
        for i in np.arange(self.nepochs):
            out = out + " {l:.3f}".format(l=self.lunation[epoch_indx[i]])
        out = out + "\n"
        out = out + " delta = "
        for i in np.arange(self.nepochs):
            out = out + " {d:.3f}".format(d=self.delta[epoch_indx[i]])
        out = out + "\n"
        out = out + " delta_min = "
        for i in np.arange(self.nepochs):
            out = out + " {s:.3f}".format(s=self.delta_min[epoch_indx[i]])
        out = out + "\n"
        out = out + " delta_max = "
        for i in np.arange(self.nepochs):
            out = out + " {s:.3f}".format(s=self.delta_max[epoch_indx[i]])
        out = out + "\n"
        out = out + " instrument = "
        for i in np.arange(self.nepochs):
            out = out + " {s}".format(s=self.instrument[epoch_indx[i]])
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
        """Define epochs based on exposure list"""
        epoch_indx = [0]
        self.nepochs = 1
        for indx in np.arange(self.nexposures - 1) + 1:
            if(self.delta[indx] > 0.):
                epoch_indx.append(indx)
                self.nepochs = self.nepochs + 1
        self.epoch_indx = np.array(epoch_indx, dtype=np.int32)
        self.epoch_nexposures = np.zeros(self.nepochs, dtype=np.int32)
        for indx in np.arange(self.nepochs - 1):
            self.epoch_nexposures[indx] = (self.epoch_indx[indx + 1] -
                                           self.epoch_indx[indx])
        self.epoch_nexposures[-1] = self.nexposures - self.epoch_indx[-1]
        return

    def evaluate_next(self, mjd_past=None, mjd_next=None,
                      lunation_next=None, check_lunation=True):
        """Evaluate next choice of observation (not well-tested)"""
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
    epoch_array(): return an ndarray with epoch-oriented list of cadences

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

        instrument : list of str
            instrument for each exposure
"""
        cadence = Cadence(*args, **kwargs)
        self.cadences[name] = cadence
        self.ncadences = self.ncadences + 1

    def check_exposures(self, one=None, two=None, indx2=None, sub1=None,
                        epoch_level=True):
        """Is exposure set in cadence two consistent with cadence one?

        Parameters:
        ----------

        one : string
            name of cadence #1

        two : string
            name of cadence #2

        indx2 : ndarray of np.int32
            exposures in cadence #2

        epoch_level : boolean
            compare sequences at epoch level not exposure level (default True)
"""
        eps = 1.e-4  # generously deal with round-off
        nexp = len(indx2)

        if(sub1 is None):
            if(epoch_level):
                sub1 = np.arange(self.cadences[one].nepochs)
            else:
                sub1 = np.arange(self.cadences[one].nexposures)

        # Check number of exposures, if at epoch level
        if(epoch_level):
            for indx in np.arange(nexp):
                if(self.cadences[one].epoch_nexposures[sub1[indx]] >
                   self.cadences[two].epoch_nexposures[indx2[indx]]):
                    return(False)

        # For the subsequent checks, convert to exposure index if we
        # are at the epoch level
        if(epoch_level):
            indx2 = self.cadences[two].epoch_indx[indx2]
            sub1 = self.cadences[one].epoch_indx[sub1]

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

    def cadence_consistency(self, one, two, return_solutions=True,
                            epoch_level=True):
        """Is cadence #1 consistent with cadence #2?

        Parameters:
        ----------

        one : string
            name of cadence #1

        two : string
            name of cadence #2

        return_solutions: boolean
            return list of solutions? (default False)

        epoch_level : boolean
            compare sequences at epoch level not exposure level (default True)

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
        if(epoch_level):
            n1 = self.cadences[one].nepochs
            n2 = self.cadences[two].nepochs
        else:
            n1 = self.cadences[one].nexposures
            n2 = self.cadences[two].nexposures

        # Check which exposures you can start on
        possibles = []
        for first in np.arange(n2 - n1 + 1):
            ok = self.check_exposures(one=one, two=two, indx2=[first],
                                      sub1=[0], epoch_level=epoch_level)
            if(ok):
                possibles.append([first])
        if(len(possibles) == 0):
            success = 0
            if(return_solutions):
                return(success, possibles)
            else:
                return(success)

        # Now find sequences starting from there
        for nsub1 in np.arange(n1 - 1) + 2:
            current_possibles = possibles
            possibles = []
            for indx in range(len(current_possibles)):
                possible = current_possibles[indx]
                remaining_start = possible[-1] + 1
                nremaining = n2 - possible[-1] - 1
                ok = 1
                if(nremaining >= n1 - len(possible)):
                    for next_possible in (remaining_start +
                                          np.arange(nremaining)):
                        try_possible = possible.copy()
                        try_possible.append(next_possible)
                        ok = self.check_exposures(one=one, two=two,
                                                  indx2=try_possible,
                                                  sub1=np.arange(nsub1),
                                                  epoch_level=epoch_level)
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

        epoch_targets : list of ndarray of np.int32
            for each epoch, list of target indices

        exposure_targets : ndarray of np.int32
            for each epoch, target index

        Notes:
        -----

        Designed to maximize total "value" of targets observed. Will
        not observe partial cadences.
"""
        ntargets = len(target_cadences)
        nepochs_field = self.cadences[field_cadence].nepochs
        nexposures_field = self.cadences[field_cadence].nexposures
        if(value is None):
            value = np.ones(ntargets)
        else:
            value = np.array(value)

        pack = []
        epochs = [0] * nepochs_field

        # Find solutions for each target
        for target_cadence in target_cadences:
            count, solns = self.cadence_consistency(target_cadence,
                                                    field_cadence,
                                                    return_solutions=True,
                                                    epoch_level=True)
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
        for target, target_cadence in zip(pack, target_cadences):
            indxs = 0
            targetvar = []
            for soln in target:
                indxe = 0
                solnvar = []
                target_indx = 0
                for epoch in soln:
                    name = "{t}-{s}-{e}".format(t=indxt, s=indxt, e=indxe)
                    if(epoch):
                        nexposures = int(self.cadences[target_cadence].epoch_nexposures[target_indx])
                        epochvar = solver.IntVar([0, nexposures], name)
                        target_indx = target_indx + 1
                    else:
                        epochvar = solver.IntVar([0], name)
                    solnvar.append(epochvar)
                    indxe = indxe + 1
                indxs = indxs + 1
                targetvar.append(solnvar)
            indxt = indxt + 1
            packvar.append(targetvar)

        # Constraint for each solution that all or none of
        # its epochs, and only those epochs, are used, and that
        # each of the epochs has exactly the right number of
        # exposures.
        for target, targetvar, target_cadence in zip(pack, packvar,
                                                     target_cadences):
            for soln, solnvar in zip(target, targetvar):
                firstvar = None
                for epoch, epochvar in zip(soln, solnvar):
                    if(epoch):
                        if(firstvar):
                            solver.Add((epochvar > 0) == (firstvar > 0))
                        else:
                            firstvar = epochvar
                    else:
                        solver.Add(epochvar == 0)

        # Constraint for each epoch that no more than
        # the total number of available exposures is taken
        for iepoch in range(nepochs_field):
            e = [solnvar[iepoch] for targetvar in packvar
                 for solnvar in targetvar]
            solver.Add(solver.Sum(e) <=
                       int(self.cadences[field_cadence].epoch_nexposures[iepoch]))

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

        # Retrieve list of targets for each epoch
        epoch_targets = [np.zeros(0, dtype=np.int32)] * nepochs_field
        if collector.SolutionCount() > 0:
            best_solution = collector.SolutionCount() - 1
            for itarget, targetvar in zip(range(ntargets), packvar):
                for solnvar in targetvar:
                    for iepoch, epochvar in zip(range(nepochs_field), solnvar):
                        if(collector.Value(best_solution, epochvar)):
                            epoch_targets[iepoch] = np.append(epoch_targets[iepoch], itarget)

        # Now convert into list of targets for each exposure
        exposure_targets = np.zeros(nexposures_field, dtype=np.int32) - 1
        iexposure = 0
        target_epoch = np.zeros(len(target_cadences), dtype=np.int32)
        for iepoch in np.arange(nepochs_field, dtype=np.int32):
            for itarget in np.arange(len(epoch_targets[iepoch]), dtype=np.int32):
                ctarget = epoch_targets[iepoch][itarget]
                ccadence = self.cadences[target_cadences[ctarget]]
                nexposures = ccadence.epoch_nexposures[target_epoch[ctarget]]
                for indx in iexposure + np.arange(nexposures):
                    exposure_targets[indx] = ctarget
                iexposure = iexposure + nexposures
                target_epoch[ctarget] = target_epoch[ctarget] + 1

        return(epoch_targets, exposure_targets)

    def fromarray(self, cadences_array=None):
        """Add cadences to ccadence list from an array

        Parameters:
        -----------

        cadences_array : ndarray
            ndarray with columns 'NEXPOSURES', 'LUNATION', 'DELTA',
            'DELTA_MIN', 'DELTA_MAX', 'CADENCE', 'INSTRUMENT'
"""
        for ccadence in cadences_array:
            nexp = ccadence['NEXPOSURES']
            instruments = [ii.decode().strip() for ii in ccadence['INSTRUMENT'][0:nexp]]
            self.add_cadence(nexposures=ccadence['NEXPOSURES'],
                             lunation=ccadence['LUNATION'][0:nexp],
                             delta=ccadence['DELTA'][0:nexp],
                             delta_min=ccadence['DELTA_MIN'][0:nexp],
                             delta_max=ccadence['DELTA_MAX'][0:nexp],
                             name=ccadence['CADENCE'].decode().strip(),
                             instrument=instruments)
        return

    def fromfits(self, filename=None):
        """Add cadences to ccadence list from a FITS file

        Parameters:
        -----------

        filename : str
            File name to read from

        Notes:
        -----

        Expects a valid FITS file with columns 'NEXPOSURES',
            'LUNATION', 'DELTA', 'DELTA_MIN', 'DELTA_MAX', 'CADENCE',
            'INSTRUMENT'
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
        cadence0 = [('CADENCE', np.dtype('a40')),
                    ('NEXPOSURES', np.int32),
                    ('DELTA', np.float64, max_nexp),
                    ('LUNATION', np.float32, max_nexp),
                    ('DELTA_MAX', np.float32, max_nexp),
                    ('DELTA_MIN', np.float32, max_nexp),
                    ('INSTRUMENT', np.dtype('a10'), max_nexp)]
        cads = np.zeros(self.ncadences, dtype=cadence0)
        names = self.cadences.keys()
        for indx, name in zip(np.arange(self.ncadences), names):
            nexp = self.cadences[name].nexposures
            cads['CADENCE'][indx] = name
            cads['NEXPOSURES'][indx] = nexp
            cads['DELTA'][indx, 0:nexp] = self.cadences[name].delta
            cads['DELTA_MIN'][indx, 0:nexp] = self.cadences[name].delta_min
            cads['DELTA_MAX'][indx, 0:nexp] = self.cadences[name].delta_max
            cads['LUNATION'][indx, 0:nexp] = self.cadences[name].lunation
            cads['INSTRUMENT'][indx, 0:nexp] = self.cadences[name].instrument
        return(cads)

    def epoch_array(self):
        """Return cadence epoches as a record array

        Returns:
        -------

        cadences : ndarray
            information on each cadence
"""
        neps = np.array([c.nepochs for c in self.cadences.values()])
        max_nep = neps.max()
        cadence0 = [('CADENCE', np.dtype('a40')),
                    ('NEPOCHS', np.int32),
                    ('NEXPOSURES', np.int32, max_nep),
                    ('DELTA', np.float64, max_nep),
                    ('LUNATION', np.float32, max_nep),
                    ('DELTA_MAX', np.float32, max_nep),
                    ('DELTA_MIN', np.float32, max_nep),
                    ('INSTRUMENT', np.dtype('a10'), max_nep)]
        cads = np.zeros(self.ncadences, dtype=cadence0)
        names = self.cadences.keys()
        for indx, name in zip(np.arange(self.ncadences), names):
            nep = self.cadences[name].nepochs
            cads['CADENCE'][indx] = name
            cads['NEPOCHS'][indx] = nep
            epoch_indx = self.cadences[name].epoch_indx
            cads['NEXPOSURES'][indx, 0:nep] = self.cadences[name].epoch_nexposures
            cads['DELTA'][indx, 0:nep] = self.cadences[name].delta[epoch_indx]
            cads['DELTA_MIN'][indx, 0:nep] = self.cadences[name].delta_min[epoch_indx]
            cads['DELTA_MAX'][indx, 0:nep] = self.cadences[name].delta_max[epoch_indx]
            cads['LUNATION'][indx, 0:nep] = self.cadences[name].lunation[epoch_indx]
            instruments = [self.cadences[name].instrument[i] for i in epoch_indx]
            cads['INSTRUMENT'][indx][0:nep] = instruments
        return(cads)
