# BabaYaga@NLO

`BabaYaga@NLO` is an event generator for Bhabha scattering, $\mu^+\mu^-$,
$\pi^+\pi^-$, photon pair production and radiative channels $\mu^+\mu^-\gamma$,
$\pi^+\pi^-\gamma$ at flavour factories currently developed by E. Budassi,
C.M. Carloni Calame, M. Ghilardi, A. Gurgone, G. Montagna, M. Moretti,
O. Nicrosini, F. Piccinini and F.P. Ucci. [(*)](#footnote)

## Compiling and running the code

Just tune the file Makefile for your needs and run `make`.
Default compilation is performed with main.F and produces the executable `babayaga`.

A template main (driver driver_gen_events.F) is provided,
containing a minimal version of a main program calling the wrapper routine to generate events.
To compile the wrapped version of the main just run

	make libbabayagafull.a
	gfortran -O3 -fPIC driver_gen_events.F libbabayagafull.a -o babayaga-wrap


## External programs

In $\textrm{BabaYaga@NLO}$, we make use of the following external programs:

- [`RECOLA`](https://recola.gitlab.io/)
- [`COLLIER`](https://collier.hepforge.org/)
- [`LoopTools`](https://feynarts.de/looptools/)
- [`RANLUX`](https://luscher.web.cern.ch/luscher/ranlux/) (C version)

## Vacuum polarization choice

In `BabaYaga@NLO`, the hadronic contribution to the vacuum polarization
can be calculated by means of three different routines:
-  `NSK` (default), Novosibirsk parameterization by Fedor Ignatov
   (http://cmd.inp.nsk.su/~ignatov/vpl/)

- `HADR5N16`, by F. Jegerlehner
   (http://www-com.physik.hu-berlin.de/~fjeger/)

- `HMNT`, by K. Hagiwara, A.D. Martin, D. Nomura and T. Teubner
   (see: Phys. Rev. D, 2004, 69:093003; Phys. Lett. B, 2007, 649:173).

You can choose it in the interactive menu. Please notice that
the `HMNT` and `HADR5N16` routines, as implemented here, are not reliable on top of the narrow resonances,
where the only reliable is `NSK`.

## Radiative channels

In the case of radiative channels $\mu^+\mu^-\gamma$ and $\pi^+\pi^-\gamma$ you can choose
between gauge-invariant subset of corrections:
- ISR corrections (Pure ISR at tree level and NLO-ISR corrections)
- FSR corrections (Pure FSR at tree level and NLO-FSR corrections)
- FULL corrections

This can be done by tuning the parameter ifisrfsr in the file `userinterface.F` with 
the following convention:
- `10` &rarr; ISR
- `01` &rarr; FSR
- `11` &rarr; FULL

Please notice that to change the subset to be considered, you have to run `make` again.

## Hadronic channels

In the case of the $\pi^+\pi^-$ channel, the form factor is introduced within three
different approaches:
- F x sQED
- GVMD
- FsQED

In the case of the  $\pi^+\pi^-\gamma$ channel, the form factor is introduced within the factorized 
approach. 
Nevertheless, both in $\pi^+\pi^-$ and $\pi^+\pi^-\gamma$, the form factor parametrizazion can be chosen
from the following options: 

-  `strong2020`
-  `cmd2`
-  `cmd3`
-  `bwsum2`
-  `bwsum3`
-  `snd`
-  `babar`
-  `besiii`
-  `kloe2`
-  `phokhara`
-  `bern`

When you run BabaYaga, enter at the prompt the variable and value you want to
change and then type `run`. You can also type `help` to print a brief
description of the input parameters and their values or `quit` to abort the
run.

## Input parameters

- `fs`		final state (ee/gg/mm/pp/pi/mr/pr)
	- `ee` &rarr; $e^+e^-\to e^+e^-$.
	- `mm` &rarr; $e^+e^-\to \mu^+\mu^-$.
	- `gg` &rarr; $e^+e^-\to \gamma\gamma$.
	- `pp` &rarr; $e^+e^-\to \pi^+\pi^-$.
	- `pi` &rarr; $e^+e^-\to \pi^+\pi^-$ (ONLY ISR CORRECTION).
	- `mr` &rarr; $e^+e^-\to \mu^+\mu^-\gamma$.
	- `pr` &rarr; $e^+e^-\to \pi^+\pi^-\gamma$.


- `ecms` &rarr; center of mass energy, in GeV

- `thmin` &rarr; minimum scattering angle for leptons (pions) in the final state, in degrees

- `thmax` &rarr; maximum scattering angle for leptons (pions) in the final state, in degrees

- `thgmin` &rarr; minimum scattering angle for the signal photon, in degrees

- `thgmax` &rarr; maximum scattering angle for for the signal photon, in degrees

- `zmax` &rarr; maximum acollinearity angle between finale state leptons/photons, in degrees

- `emin` &rarr; minimum energy for leptons in the final state, in GeV

- `egmin` &rarr; minimum energy for the signal photon, in GeV

- `massmin` &rarr; minimum invariant mass of the $l^+l^-$ system, in GeV

- `massmax` &rarr; maximum invariant mass of the $l^+l^-$ system, in GeV

- `nev` &rarr; number of events to be generated

- `path` &rarr; directory where to store output files. In *nix systems, the directory is automatically created.

- `saveevents` &rarr; if saving an ascii file where weighted or unweighted (according to `mode`) events are save. The file is `path/events.dat`.

- `iffpi` &rarr; the way $F_\pi(q^2)$ is introduced in the calculation:
	- For radiative $\pi^+\pi^-\gamma$
 		- `0` &rarr; Pion form factor off
		- `1` &rarr; Pion form factor $\text{F}\times\text{sQED}$

	- For $\pi^+\pi^-$ production
		- `0` &rarr; Pion form factor off
		- `1` &rarr; Pion form factor $\text{F}\times\text{sQED}$
		- `2` &rarr; Pion form factor GVMD
		- `3` &rarr; Pion form factor FsQED

- `what_ffpi` &rarr; to set which parametrization of the pion form factor must be used.

- `arun`	  &rarr; sets $\alpha(s)$ routine:
 	- `off`   &rarr; sets alpha running off
	- `nsk`   &rarr; `NSK` parameterization by Fedor Ignatov
	- `hadr5` &rarr; `HADR5N16` routine
	- `hmnt`  &rarr; `HMNT` routine


- `mode` &rarr; sets if the requested number of events (nev) are weighted or unweighted

- `eps` &rarr; sets the soft/hard photon energy separator, in ecms/2 units.
		Results are completely independent from its choice (provided it is small)

- `ord` &rarr; it sets which photonic radiative corrections are included:
	- `born` means Born cross section,
	- `alpha` means cross section at order alpha
	- `exp` is the best, exponentiated cross section
	 
- `model` &rarr; it sets the model for radiative corrections.
	- `matched` is the best one
	- `ps` is very similar to older `BabaYaga@NLO` releases, also in this case the ps is active only for non radiative process.
	 
- `seed` &rarr; an integer number to initialize the random number generator

- `nphot` &rarr; only a fixed number of `nphot` (hard) photons are generated. A negative value means all possible photons

- `nwrite` &rarr; output files in `path` are written every `nwrite` events, if `nwrite` is set to a negative value, output is written approximately every `-nwrite` seconds

- `nsearch` &rarr;	nsearch events are generated to find the maximum value of the cross section, after which also events unweightening is started.

- `verbose` &rarr; it toggles some verbose output, only for debugging

- `sdmax` &rarr; the starting maximum value for the cross section

## User modifiable routines
 
The subroutines the user may need to modify are:
- `cuts` (in the file `cuts.F`), used to impose event selection
	criteria on the generated events. Cuts are always performed
	in the lab reference frame.

- `initstorage` (in the file `routines.F`): it initializes the environment
	to store unweighted/weighted events in the `events.dat` file.

- `eventstorage` (in the file `storage.F`): it stores the
	unweighted/weighted events in the event file. By default, it saves the
	four-moment of the final state $l^+,l^-$, and the emitted photons.
	Any other event variable can be of course saved.

## Output files

The output files are saved in the `path` directory. The files are
- `events.dat`: it is the file where unweighted/weighted events are stored, if `saveevents` is set to `yes`.

- `statistics.txt`: it is the file (dumped every `nwrite` points) where cross sections, statistics information,
  input parameters, etc. are printed. After the input parameters entered by the user, the weighted
	integrated cross section is printed. It is subdivided by photon multiplicity: the single cross sections
	depend on the `eps`  parameter but the sum does not. After weighted cross section information, the
	unweighted events statistics is reported. The user must pay attention to the "biases" cross sections:
	they account for the biases due to negative weights events (highly improbable, we never saw them) and
	unweightening failure due to an under-estimated maximum of the differential cross section. The
	latter bias is usually negligible, its size is estimated by the corresponding cross section.

-	distribution files: in these files some differential distribution,
	calculated with weighted events, are written. The file names should be self-explaining and the data
	are in the form (e.g. to be plot easily with `gnuplot`)
	lower bin edge - differential cross section - corresponding error
  
 
## Input card and default cuts
An input card (`input_rad`) that generates 10k events (with arbitrary cuts) at order
exp for the radiative process $\mu^+\mu^-\gamma$ is provided.

The authors

Andrea, Carlo, Ettore, Francesco, Fulvio, Guido, Marco, Mauro, Oreste
<br>
<br>

#### Footnote
Former developers include G. Balossini, L. Barze`and C. Bignamini.
