# dynamics-w90: Time-dependent dynamics and band properties from Wannier functions #

This is a small collection of programs for:

* computing time-depndent dynamics in the presence of electric fields 
* computing band properties including orbital weights, spin texture, Berry curvature and orbital angular momentum.

The programs are written in modern Fortran, while interface and plotting scripts are written in python.

## Directory structure ##

The repository includes the following directories:

### inp ###

Input files for the programs go here. We use the Fortran namelist format for most input.
We have included some example input files for monolayer WSe2.

### data ###

Here we have included output files from `Wannier90` that are read by the programs.

### examples ###

Some example scripts for plotting the calculated band properties.

### exe ###

The executables will be placed here after successfull compilation.

### python_utils ###

Some useful python modules for reading output of the programs and for plotting.

### out ### 

Output files in hdf5 format will be stored here.

## Installation ##

### Dependencies ###
* [CMake](https://cmake.org)
* [HDF5](https://www.hdfgroup.org). This is optional if HDF5 support is requested.
* A version of BLAS/LAPACK
* For the MPI version of the time-evolution code (`wann_evol_mpi.x`): `mpi`. This is optional. 
* For the interface and plotting scripts, we use `numpy`, `matplotlib`. For reading `hdf5` output, `h5py` is required.

### Compilation ###

For installing the `dynamics-w90` programs, create a configure script or execute in the terminal:

    # From the root directory .../dynamics-w90/
    # Create build directory
    mkdir build
    cd build

    # Cmake configuration step (replace paths according to your system)
    FC=[Fortran compiler] \
    cmake \
        -DCMAKE_INSTALL_PREFIX = install_path \
        -DCMAKE_BUILD_TYPE=[Release|Debug] \
        -DCMAKE_INCLUDE_PATH=/opt/local/include \
        -DCMAKE_LIBRARY_PATH=/opt/local/lib \
        -DCMAKE_Fortran_FLAGS="[optimization flags|debug flags]" \
        ..

    # compile
    make

After commpleted compilation, the `exe` directory will be created and the executables will be placed there. 

HDF5
====

To enable HDF5 support use the `hdf5` option in cmake

    -Dhdf5=ON

Without HDF5 support the programs will write output to txt files.

BLAS/LAPACK
===========

Specify your version of BLAS/LAPACK by 

    -Dlapackblas_libraries="lib"

For instance, to link against [OpenBlas](https://github.com/xianyi/OpenBLAS), use

    -Dlapackblas_libraries="-lopenblas"

OpenMP
======

To enable OpenMP support use the `omp` option in cmake

    -Domp=ON

MPI
===

To enable MPI support use the `mpi` option in cmake

    -Dmpi=ON

In this case the Fortran compiler must be an MPI compiler like `mpif90`, `mpifort` or `mpiifort`.

## Running the programs ##

### wann_evol.x ###

Run 

    ./exe/wann_evol.x inp/input.inp out/prefix

The input file is a Fortran namelist file with the following input variables:

#### SYSPARAMS ####

* `MuChem`: the chemical potential (a.u.)
* `FixMuChem`: if `.true.`, the input chemical potential will be used. Otherwise the chemical potential will be recalculated to match the given filling.
* `Filling`: Filling of the bands, corresponding to the total number of electrons per unit cell (per spin without SOC).
* `Beta`: inverse temperature (a.u.)
* `Nk1,Nk2,Nk3`: Discretization of the Brillouin zone in the three directions. Set `Nk3=1` for 2D systems.
* `file_ham`: File with the Wannier Hamiltonian. This is the `_tb.dat` file obtained from `Wannier90`.
* `gauge`: Velocity gauge (`gauge=0`), dipole gauge (`gauge=1`), empirical velocity gauge (`gauge=3`), Peierls substitution (`gauge=4`)

#### TIMEPARAMS ####

* `Nt`: The number of time steps.
* `Tmax`: The propagation time (a.u.)
* `output_step`: Calculate observables every `output_step` time steps.

#### FIELDPARAMS ####

* `ApplyField`: If `.true.`, the laser field will be read from file. Otherwise the field-free evolution will be computed.
* `file_field`: Data file with three-dimensional electric field. The following format is expected: time (1st column), Ex (2nd column), Ey (3rd column), Ez (4th column). The time and field strength is expected in atomic units.

If compiled with HDF5 support, after running the code the `hdf5` file `out/prefix_observables.h5` will be produced. The observables can be plotted by using the script `python_utils/plot_obs.h5`.
Without HDF5 support there will be several output files `out/prefix_etot.txt`, `out/prefix_curr.txt` etc. that can directly be plotted or read by numpy's `loadtxt`.

### wann_evol_mpi.x ###

The input variables and output are identical to `wann_evol.x`. 

### wann_calc.x ###

This program computes band properties. Run 

    ./exe/wann_calc.x inp/input.inp out/prefix

The input file is a Fortran namelist file with the following input variables:

#### CALCOPT ####

This namelist controls which properties will be computed and store to file. 

* `calc_orbweight`: Triggers the calculation of the orbital weight of each Wannier function.
* `calc_spin`: Triggers the calculation of the three-dimensional spin texture. Works only if the Hamiltonian is in spin space and the spin quantization axis is in z direction.
* `calc_berry`: Triggers the calculation of the Berry curvature.
* `calc_oam`: Triggers the calculation of orbital angular momentum (OAM).
* `calc_evecs`: If `.true.`, the complex eigenvectors will be written to file.
* `gauge`: Velocity gauge (`gauge=0`) or dipole gauge (`gauge=1`) for the calculation of Berry curvature or OAM.

#### HAMILTONIAN ####

* `file_ham`: File with the Wannier Hamiltonian. This is the `_tb.dat` file obtained from `Wannier90`.
* `w90_with_soc`: Tells the program that the Hamiltonian has been obtained with SOC.

#### KPOINTS ####

* `kpoints_type`: Controls the format of the k-points input. Options are 1) `"path"`, 2) `"list"`, or 3) `"grid"`
* `file_kpts`: The input file specifying the k-points. Used if `kpoints_type="path"` or `kpoints_type="list"`. 
* `Nk1,Nk2,Nk3`: Discretization of the Brillouin zone in the three directions. Only used if `kpoints_type="grid"`.

If `kpoints_type="list"` the code expects a list of k-points with kx (1st column), ky (2nd column), and kz (3d column). Input in reduced coordinates is expected.

If `kpoints_type="path"` a path in k-space will be constructed from the input. The `file_kpts` file has the following format:

    npoints ndim
    nseg1 nseg2 ...
    point1 
    point2
    ..

Here, `npoints` is the number of points to pass through, while `nseg1` is the number of segments between `point1` and `point2` and so on. Below the special points are listed.

## Plotting ##

Check the scripts in the `examples` directory. Scripts ending with `_txt.py` do not require HDF5 support and directly read txt-based output.




