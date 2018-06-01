A basic task that needs to be performed is to decide what cadence to
take in each field. This task is closely related to deciding which
fibers to assign to which target for each design. 

Here we treat the cadences as defined by a set of exposures, each of
which has a "day" and a "softness" associated with it. The first
exposure on the list has day=0 (though in principle, other exposures
can have negative "day" values).

Cadence Consistency Matrix
==========================

For a given field, it will have some cadence. For each robot in the
field, it can reach certain targets, which have some cadence
requirements.

We need to know the answer to the question:

 "Is cadence i achievable in a field with cadence j?"

where we might as well find out the answer for all pairs of cadences
under consideration.

For example, a single exposure cadence is always achievable, but if we
need two exposures on a target in cadence i then its achievability
depends on the requirements on their cadence and what is on offer
under cadence j. That is, if cadence i needs a month separation, but
cadence j has yearly separations, then it won't work. 

If we can determine this, we can construct a "cadence consistency
matrix", C_{ij} which contains the answer to the above question
encoded as a 0 or 1 (i.e. False or True). Note that this matrix is
decidedly not symmetric. This matrix will be the basis on which we
make decisions about which targets the robots in a given field with a
chosen cadence should look for.

To determine the cadence consistency matrix, we have to look at the
individual exposures. Imagine we want to get exposures in the cadence
specified in cadence i. Well, the first exposure in the list needs to
correspond to some exposure in cadence j. So we will check each
possibility in turn.  

For each choice of day_i(0) = day_j(jstart), we now create an Exposure
Consistency Matrix. The exposure consistency matrix contains the
answer to the question for each i', j':

 "Is exposure i' consistent with exposure j' offset by day_j(jstart)?"

The Exposure Consistency Matrix, E_{ij} is 0s and 1s. The condition
that the cadence i is achievable under cadence j is equivalent to
being able to construct a matrix w_{ij} for which the following is
true for some choice of jstart:

 w_{ij} <= E_{ij} for all i, j, AND 
 \sum_i w_{ij} <= 1 for all j, AND
 \sum_j w_{j} == 1 for all i

The w_{ij} are the "solutions" -- the choices of which exposures for
cadence i to put in which exposures for cadence j.

The solutions to the Exposure Consistency Matrix problem can be found
with standard techniques in constraint programming. Thus, we can
determine the Cadence Consistency Matrix.

Packing Targets into a Cadence
==============================

For each robot in a field observed under a certain cadence, there are
then some set of targets available to it. We want then to pack those
targets into the field cadence in such a way as to maximize the
science (i.e. by getting the most targets, or something). 

We can solve this problem in a similar way to solving the exposure
consistency matrix. In this case we construct a matrix A_{ijk}, where
i indexes the target, j indexes all the patterns that solve the
the ECM for the target i cadence under the field cadence, and k
indexes the epoch. Note this isn't quite a 3D matrix, because there
are different numbers J of patterns for each target cadence. A_{ijk}
is tells you whether for target i, under pattern j, at epoch k, the
target would be observed. 

Then you need to find w_{ijk} such that:

 w_{ijk} <= A_{ijk} 
 \sum_{ij} w_{ijk} <= 1 (the robot can only observe one target each epoch)
 \sum_{j} [(\sum_k w_{ijk}) > 0] <= 1 (only one cadence j chosen for a target)
 w_{ijk} = w_{ijk'} for any choices k, k' (all epochs in a cadence, or none)

If you assign a value to each target i, you can define a total value:

 V = \sum_i v_i [(\sum_jk w_{ijk}) > 0] 

Constraint programming techniques can then find solution w_{ijk} that
maximize the value, and you then known for each epoch k which target i
to observe. 
