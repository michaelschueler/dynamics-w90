project: dynamics-w90 
src_dir: ./src
	./scitools
output_dir: ./doc
exclude: lebedev_weights.F90 
	units_inc.f90 
	formats.h
exclude_dir: ./scitools/math/lebedev/
author: Michael Schueler
author_description: The `dynamics-w90` package: Berry phase properties, light-induced dynamics, photoemission & more from Wannier functions.
project_github: https://github.com/michaelschueler/dynamics-w90
project_download: https://github.com/michaelschueler/dynamics-w90/releases/latest
github: https://github.com/michaelschueler
display: public
         protected
         private

Description
-----------------

**dynamocs-w90** is a software package for computing Berry phase properties, light-induced dynamics, photoemission & more from Wannier functions.

Here is a list of the main programs:

1. [[wann_evol]] / [[wann_evol_mpi]]: Computes the time-dependent dynamics of an electron system described by a Wannier Hamiltonian upon laser excitation and computes various observables. Direct implementation of the gauge-invariant formulation from [Phys. Rev. B 103, 1155409 (2021)](https://link.aps.org/doi/10.1103/PhysRevB.103.155409) 
2. [[wann_calc]]: Calculation of band properties including orbital weights, spin texture, Berry curvature and orbital angular momentum (OAM). We implemented the modern theory of OAM to include all non-local contributions.
3. [[wann_prune]]: Compresses a Wannier Hamiltonian by cutting small hopping amplitudes. Output can be written in original Wannier90 or in hdf5 format.
4. [[wann_soc]]: Extends a Wannier Hamiltonian in spin space and adds atomic spin-orbit coupling (SOC).
5. [[arpes]] / [[arpes_mpi]]: Computes angle-resolved photoemission spectrum from Wannier functions. 

@Note The code is under active development. Other functionalities will be added soon!

Official Releases
-----------------

The **current stable release** can be [downloaded
on GitHub](https://github.com/michaelschueler/scitools/releases/latest). The
documentation for the current version is included in the the release.

A list of all past releases, links to their documentation, and the
change log can be found on the
[releases page](https://github.com/michaelschueler/scitools/releases/index.html).
