# Packmol

Packmol - Creates Initial Configurations for Molecular Dynamics Simulations

This page contains the version history of Packmol. You can download packmol from this page, but give preference to the official Packmol page:

http://m3g.iqm.unicamp.br/packmol

Detailed installation and usage instructions are available at this home-page.

## What is Packmol

Packmol creates an initial point for molecular dynamics simulations by packing molecules in defined regions of space. The packing guarantees that short range repulsive interactions do not disrupt the simulations.

The great variety of types of spatial constraints that can be attributed to the molecules, or atoms within the molecules, makes it easy to create ordered systems, such as lamellar, spherical or tubular lipid layers.

The user must provide only the coordinates of one molecule of each type, the number of molecules of each type and the spatial constraints that each type of molecule must satisfy.

The package is compatible with input files of PDB, TINKER, XYZ and MOLDY formats.

See http://m3g.iqm.unicamp.br/packmol for more information.

## Installation instructions

### From source:

1. Download the `.tar.gz` or `.zip` files of the latest version from: https://github.com/m3g/packmol/releases

2. Unpack the files, for example with: 
   ```
   tar -xzvf packmol-20.12.0.tar.gz
   ```
   or
   ```
   unzip -xzvf packmol-20.12.0.zip
   ```
   substituting the `20.12.0` with the correct version number.

3. Go into the directory, and compile the package (we assume `gfortran`) is installed in your system:
    ```
    cd packmol
    make
    ```

4. An executable called `packmol` will be created in the main directory. Add that directory to your path.

### Using the Fortran Package Manager

packmol now is compatible with the Fortran Package Manager `fpm`, instructions
of installing `fpm` can be found in:
https://fpm.fortran-lang.org/en/install/index.html#install

1. Download the `.tar.gz` or `.zip` files of the latest version from: https://github.com/m3g/packmol/releases

2. Unpack the files, for example with: 
   ```
   tar -xzvf packmol-20.12.0.tar.gz
   ```
   or
   ```
   unzip -xzvf packmol-20.12.0.zip
   ```
   substituting the `20.12.0` with the correct version number.

3. Go into the directory, to install the package with `fpm` just run:
   ```
   fpm install --profile release
   ```
   this will compile and send the executable somewhere in your `PATH`.
   By default (on Linux systems) it will be `~/.local/bin`. Making it available
   as a `packmol` command anywhere in your computer.

   `fpm` will look for Fortran compilers automatically and will use `gfortran`
   as default. To use another compiler modify the environment variable
   `FPM_FC=compiler`, for example for `ifort` `FPM_FC=ifort`.

## Usage

An user guide, examples, and tutorials, are available at: http://m3g.iqm.unicamp.br/packmol

## References

Please always cite one of the following references in publications for which Packmol was useful:

L Martinez, R Andrade, EG Birgin, JM Martinez, Packmol: A package for building initial configurations for molecular dynamics simulations. Journal of Computational Chemistry, 30, 2157-2164, 2009. (http://www3.interscience.wiley.com/journal/122210103/abstract)

JM Martinez, L Martinez, Packing optimization for the automated generation of complex system's initial configurations for molecular dynamics and docking. Journal of Computational Chemistry, 24, 819-825, 2003.
(http://www3.interscience.wiley.com/journal/104086246/abstract)



