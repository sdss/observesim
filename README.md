# Observing Simulation Code for SDSS-V

This product is for observing simulation code for SDSS-V. This
location on GitHub is temporary before being transferred to the main
SDSS repositories.

## Master Schedule

The `master` module allows the interpretation of the master schedule and
the calculation of basic astronomical parameters for each night within
it.

The master schedule itself is kept as a Yanny file at:
```
 $OBSERVESIM_DIR/data/master_schedule.par
```
and an example is given below.  It contains a list of events labeled
by local time and date, each of which either turns the survey on or
turns it off.

The list of observatories is at 
```
 $OBSERVESIM_DIR/data/observatories.par
```
and an example is given below.

An example of how to use the master schedule module is as follows:

```
import observesim.master as master

schedule = master.Master()
template = "mjd={mjd}, illumination={illumination}"
for mjd in schedule.mjds:
    night = master.Night(mjd=mjd)
    illumination = night.moon.illumination(mjd=mjd)
    print(template.format(mjd=mjd, illumination=illumination))
```

The `Night` class also has attributes which are the floating point
MJDs of evening and morning twilight.  Thus, a simulation can use this
module to loop through all of the nights, and then simulate each
night. The `sun` and `moon` attributes of `Night`, and the `lst()`
method, will give appropriate information for planning the night.

Dependencies:
```
 numpy
 scipy
 astropy
 pydl
 PyAstronomy
```
The reason for PyAstronomy is that it has Meeus-based routines for
astronomical calculations that are substantially faster than the more
accurate routines found in astropy. There are a couple of problems in
it with the Sun position calculation for large vectors. This part of
our code should be cleaned up

### Example contents for master_schedule.par

```
# Timezone offset in hours to apply to get to TAI
# (i.e. Greenwich time)
to_tai 7  # Mountain Standard Time

# Whether events start ("on") or stop ("off") observing
START_SURVEY on
END_SURVEY off
START_SHUTDOWN off
END_SHUTDOWN on

typedef enum {
  START,
  END,
  START_SHUTDOWN,
  END_SHUTDOWN
} EVENT;

typedef struct {
  char date[10];
  char time[5];
  EVENT event;
} SCHEDULE;

SCHEDULE 2020-07-01 12:00 START_SURVEY
SCHEDULE 2020-07-10 12:00 START_SHUTDOWN
SCHEDULE 2020-08-20 12:00 END_SHUTDOWN
SCHEDULE 2021-07-01 12:00 END_SURVEY
```

### Example contents for observatories.par

```
typedef struct {
  char observatory[100];
  double longitude;
  double latitude;
} OBSERVATORY;

OBSERVATORY APO -105.82027778 32.7797556
```
