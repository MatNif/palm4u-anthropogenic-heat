
# The kpp4palm preprocessor

---

## Introduction

`kpp4palm` is a preprocessor that creates the file `chem_gasphase_mod.f90`, where the gas phase chemistry rate equations are solved within PALM.


`kpp4palm` is based on the original unchanged Kinetic PreProcessor KPP (Damian et al., 2002, Sandu et al., 2006), [Release 2.2.3 from November 2012](http://people.cs.vt.edu/~asandu/Software/Kpp/Download/kpp-2.2.3_Nov.2012.tar.gz) and an adapted version of the KPP postprocessor KP4 (Jöckel et al, 2010), which converts the KPP-generated code to a subroutine for PALM. The adapted version of KP4 is named kpp4palm.

KPP creates code for a box model from a list of chemical reactions, which must be written in a format that can be processed by KPP. This code is converted to a module for PALM by kpp4palm.

Besides the standard scalar version of the code, also a vectorized version of chem_gasphase_mod.f90 can be generated. However, only the different flavors of the Rosenbrock solvers have been vectorized, all other KPP solvers have to run in scalar mode.

The first version of this interface to PALM (still named kp4 at that time) was created by Klaus Ketelsen in November 2016 on the basis of his previous development of the handling of KPP in MESSy2 as described by Jöckel et al. (2010).

## Directory structure

The kpp4palm preprocessor is located in directory `packages/chemistry/kpp4palm`. 

Contents of `packages/chemistry/kpp4palm` are

- Directory `docs`: Documentation (this file is in subdirectory `docs/content/Guide`).
- Directory `kpp`: KPP preprocessor creating code (Fortran in our case) from a list of chemical reactions.
- Directory `mechanisms`: Contains sub-directories with the input for KPP for some sample mechanisms and the already processed `chem_gasphase_mod.f90`).
- Directory `scripts` containing `kpp4palm.sh` for running kpp4palm (linked to `build/bin/kpp4palm` by the PALM installer).
- Directory `src` containing the kpp4palm code.
- Directory `templates` contain templates which are included into `chem_gasphase_mod.f90` when kpp4palm is run.
- Optionally, the directory `tmp_kpp4palm` can be created when running `kpp4palm`. This directory contains intermediate files, e.g. the original KPP input and output.

Output of `kpp4palm` is the file `chem_gasphase_mod.f90` which contains the Fortran code of the chemistry subroutines for a user defined chemical mechanism. This file is saved in the directory the user chooses with the option `-o`.

The required input files of `kpp4palm` and KPP are located in the directories `mechanisms/def_<mech>`, where `<mech>` stands for the name of any mechanism. A few sample mechanisms are already supplied in `mechanisms/`. More mechanisms may be added and can also be added by the user.

Each of the `def_<mech>` directories contains the following KPP input files:

- `chem_gasphase_mod.kpp` contains some instructions for KPP, such as the output language, directives for the photolysis reactions and values of the compounds which are referred as 'fixed species'. Fixed species are compounds which are usually abundant and do not vary with time on the scale of tropospheric chemistry, e.g. O2 or N2. For some mechanisms also CO2 or methane are considered as fixed, i.e. which  compounds considered as 'fixed' depend also on the mechanism. Water vapor (H2O) is always considered as a 'fixed species' in the chemistry routines, since its concentration is calculated in the meteorological part of PALM-4U (if it were not considered as 'fixed', it would be transported twice).
- `<mech>.spc` containing a list of the chemical species in KPP notation.
- `<mech>.eqn` containing a list of the chemical reactions in KPP notation.
- `UserRateLaws.f90` contain the rate laws which are actually used (currently `UserRateLaws.f90` are identical for all mechanisms) . This file is a copy of [kpp/util/UserRateLaws.f90](https://gitlab.palm-model.org/releases/palm_model_system/-/blob/master/packages/chemistry/kpp4palm/kpp/util/UserRateLaws.f90_orig) (`UserRateLaws.f90` is one of the ‘auxiliary files’ which are mentioned in the KPP documentation).
- In addition, each `mechanisms/def_<mech>` contains an already set output file `chem_gasphase_mod.f90` just in case that KPP cannot be run on a user's system (usually due to missing requirements). However, only preprocessed files for the scalar mode are supplied here.

Sample files for the vector mode are not supplied since the optimum vector length depends on the computer which is used.

The source code of PALM also contains already a file named `chem_gasphase_mod.f90`. By default, this is the mechanism phstatp, i.e. photostationary equiöibrium between ozone, NO, and NO2 plus one passive tracer named MP10. If you are not sure which mechanism is used there: The third line of `chem_gasphase_mod.f90` indicates the mechanism.

If someone wants to switch to another mechanism than the one which is included in the source code of PALM either run `kpp4palm` as described below or copy the already prepared `chem_gasphase_mod.f90` from the respective `def_<mech>` directory into USER_CODE of your JOBS/run_descriptor.

## Requirements

- FLEX library
- BISON parser generator

## Installation

- Make sure that all requirements are installed.
- Installation of KPP and `kpp4palm` is included in the standard installation of PALM (i.e. by executing the command `bash install -p <install_prefix>`).

## How to generate code for available mechanisms

You can run `kpp4palm` as follows:
```bash
kpp4palm [-h] -m <mechanism_name> -o <path> [-i <method_number>] [-v] [-l <vector_length>] [-k] [-u]
```
- Option `-h` show a help message.
- Option `-m <mechanism_name>` permits the choice of the chemical mechanism. If it is not specified, the default mechanism `phstatp` will be used.
- Option `-o <path>` sets the output directory for the generated code. This option is mandatory.
- Option `-i <method_number>` is optionally set to (method_number=0,1,2) and optimizes a part of the code by replacing indirect addressing of arrays by a sequence of statements without indirect addresses as described by Jöckel et al., 2010. If n is set to 0, then the code is not optimized. Default is 2.
- Option `-v` switches on the generation of the vector version of `chem_gasphase_mod.f90`. A vector length must be specified by using option `-l`.
- Option `-l <vector_length>` can only be applied in combination with option `-v`. See above.
- Option `-k` stands for "keep" and determines whether the temporary working directory `tmp_kpp4palm` is kept or deleted after termination of `kpp4palm`. The directory `tmp_kpp4palm` is removed when this option is omitted.
- Option `-u` stands for "update" and determines whether the output file will also be copied into the `def_<mech>` directory. This option should be applied with caution since the original file `chem_gasphase_mod.f90` in `def_<mech>` will be overwritten. The default setting is off.

During runtime a temporary directory `tmp_kpp4palm` is created. The newly created output file is copied from the temporary working directory `tmp_kpp4palm` to the output directory provided by option `-o`. If there is an already existing file `chem_gasphase_mod.f90` in the output directory, it is moved to chem_gasphase_mod.f90.sav and will be overwritten by the next run of `kpp4palm`. If the `-u` option is applied, the output file will also be copied into the `mechanisms/def_<mech>` directory.

## How to apply kpp4palm for a new mechanism

If you are not familiar with KPP, read the [KPP documentation](http://www.cs.vt.edu/~asandu/Software/Kpp/Download/kpp-2.1_UsersManual.pdf)(also locally available in `kpp4palm/kpp/doc/kpp_UserManual.pdf`) and have a look into the files of the already existing `def_<mech>` directories. The following steps are needed to create new mechanisms that you can name as you like. Choose a name and replace every occurance of <mech> in the following steps with that name:

1. Create an new subdirectory `def_<mech>` in directory mechanisms.
2. Put your new `<mech>.spc` and `<mech>.eqn` into that new directory. Photolysis frequencies must be named according to the following examples: `phot(j_no2)`, `phot(j_hcho)`, `phot(j_o3)`.
3. Copy a `chem_gasphase_mod.kpp` file into the directory and adapt it:
   - Adapt the name of the mechanism in the two `#include` statement
   - Adapt in `#INLINE F90_DATA`  the number of photolysis frequencies `nphot`
   - Adapt/extend  in `#INLINE F90_DATA` the indices in the `INTEGER, PARAMETER,PUBLIC` statement
   - Adapt/extend in #INLINE F90_DATA the character array phot_names: Note that the order of phot_names and the indices must match. Please note that the names are case sensitive. The available photolysis frequencies can be found in `chem_photolysis_mod.f90` (array `names_s`).
   - Adapt in the `#INLINE F90_INIT` section the 'fixed' species exactly to number of compound which are required for  your mechanism. Please note that water vapor is considered as fixed within the chemistry module, as it is computed somewhere else.

## Documentation and References

A local copy of the official [KPP documentation](http://www.cs.vt.edu/~asandu/Software/Kpp/Download/kpp-2.1_UsersManual.pdf) is found in `kpp4palm/kpp/doc/kpp_UserManual.pdf`.


- KPP web page: "[http://people.cs.vt.edu/asandu/Software/Kpp/](http://people.cs.vt.edu/asandu/Software/Kpp/)"
- Damian, v., A. Sandu, M. Damian, F. Potra, and G.R. Carmichael: ``The Kinetic PreProcessor KPP -- A Software Environment for Solving Chemical Kinetics'', Computers and Chemical Engineering, Vol. 26, No. 11, p. 1567-1579, 2002.
- Jöckel, P., Kerkweg, A., Pozzer, A., Sander, R., Tost, H., Riede, H., Baumgaertner, A., Gromov, S., and Kern, B.: Development cycle 2 of the Modular Earth Submodel System (MESSy2), Geosci. Model Dev., 3, 717-752, https://doi.org/10.5194/gmd-3-717-2010, 2010.
- Sandu A., and R. Sander. "[Technical Note: Simulating chemical systems in Fortran90 and Matlab with the kinetic preprocessor KPP-2.1](http://www.atmos-chem-phys.org/acp/6/187/)", Atmospheric Chemistry and Physics, Vol. 6, p. 187-195, (2006).

---

---

# Additional notes (background information about KPP and kpp4palm)

---

## General remark

KPP and kpp4palm are strongly case sensitive. Adaptations to coding conventions for PALM should therefore only be applied after the essential processing is finalized.


## How to add or modify expressions for rates

The respective files `mechanisms/def_<mech>/UserRateLaws.f90` contain the rate laws which are actually used (currently `UserRateLaws.f90` are identical for all mechanisms) . This file is a copy of [kpp/util/UserRateLaws.f90](https://gitlab.palm-model.org/releases/palm_model_system/-/blob/master/packages/chemistry/kpp4palm/kpp/util/UserRateLaws.f90_orig) (which is one of the ‘auxiliary files’ mentioned in the KPP documentation).


Further rate laws may be added in `mechanisms/def_<mech>/UserRateLaws.f90`, e.g.
```fortran
REAL(kind=dp) FUNCTION ARR2( A0,B0, TEMP )
   REAL(kind=dp) :: TEMP
   REAL(kind=dp) A0,B0
   ARR2 = A0 * EXP( -B0 /TEMP )
END FUNCTION ARR2
```

It can be extended by further rate laws. When `kpp4palm` is run, `mechanisms/def_<mech>/UserRateLaws.f90` is copied to `kpp/util/UserRateLaws.f90`.

In order to make `kpp4palm` include `ARR2` into `chem_gasphase_mod.f90`, `kpp4palm/scripts/kpp4palm.sh` is parsing `UserRateLaws.f90` for the expression `FUNCTION` and adds the names of all functions in `UserRateLaws.f90` to the list of subroutines to be processed in the file `KPP_SUBROUTINE_LIST`.

**Important:** Within `UserRateLaws.f90` the effective code lines of the subroutine must be the after any comment lines (otherwise kpp may create a memory fault).


## About the type of rate ‘constants’

Rate constants in `<mech>.eqn` (<mech> stand for any name of a mechanism) are handled differently, depending whether they are simple numbers or not.

In `kpp/src/scanner.c` there are three types of rates distinguished in `StoreEquationRate` (lines 578 ff): `NUMBER`, `EXPRESION`, and `PHOTO`.

`NUMBER` is clear. If rate is a number, then it is put to the initialization as a ‘constant rate coefficient’.
If the rate in `<mech>.eqn` includes anything which is different from a number (like brackets, `_dp`, or anything else), then the rate is of type `EXPRESION`. This type of rates is put into `UpdateRconst`.

`PHOTO` is identified by the occurrence of ‘hv’ in the reaction rate equations in `<mech>.eqn`:
`if(EqNoCase(spname,"HV")) isPhoto = 1;` within `scanner.c` (line 692).

Appending `_dp` (which is requested by the PALM team) at the end of each number makes a number to an expression and the rate is put into `SUBROUTINE UpdateRconst` (which is what we want anyway).


## The Update_Rconst calls issue

As the Box version of KPP generated code can also run for several hours, an updating of the rate constants is necessary. However, this is not required, when just a time step of a dynamical model must be covered. Then it is only necessary to call `Update_Rconst` only once at the beginning of each time step. Furthermore, `Update_SUN` is not necessary as photolysis will be calculated outside of the chemistry module.

The call of `Update_Rconst` and `Update_SUN`is removed from the code from KPP by the following lines in `kpp4palm/src/fortran_file.C` :
```c
//  Update_RCONST has only to be called once per outer timeloop in KPP_FOR_PALM

    if(ip->get_token(0) == "CALL" && ip->get_token(1) == "Update_RCONST" ) {
      lo_line.insert(0,"!DELETE ");
    cout << lo_line << endl;
    }

//  Update_SUN must not be called within in KPP_FOR_PALM

    if(ip->get_token(0) == "CALL" && ip->get_token(1) == "Update_SUN" ) {
      lo_line.insert(0,"!DELETE ");
    cout << lo_line << endl;
    }
```

## Modifications for photolysis

In file `kpp4palm/src/create_kpp_module.C`, some additional lines were added after:
```c
void create_kpp_module::create_kpp_integrate()
```

Specification of indices for the photolyis frequencies and their names must be given in `mechanisms/def_<mech>/chem_gasphase_mod.kpp`. Names of the photolysis frequencies must match with the names of the available photolysis frequencies in `chem_photolysis.f90` in the PALM source code. The following example shows the setting for the 'smog' mechanism:

```fortran
#INLINE F90_DATA
  !   Declaration of global variables for photolysis from INLINE
  INTEGER, PARAMETER :: nphot = 2
  !   phot Photolysis frequencies
  REAL(kind=dp) :: phot(nphot)

  INTEGER, PARAMETER,PUBLIC :: j_no2 = 1
  INTEGER, PARAMETER,PUBLIC :: j_rcho = 2

  CHARACTER(LEN=15), PARAMETER, DIMENSION(NPHOT) :: PHOT_NAMES =  (/ &
     'J_NO2          ','J_RCHO         '/)
#ENDINLINE
```
The declaration `REAL(kind=dp) :: phot(NPHOT)` does not really fit here, as this does not depend on the mechanism, but so far we did not find a better place.


## Modifications for fixed species

Fixed species (when necessary) are initialized in the `#INLINE F90_INIT` section of `mechanisms/def_<mech>/chem_gasphase_mod.kpp` as follows (example for 'smog' mechanism):

```fortran
#INLINE F90_INIT
  fix(indf_h2o) = qvap
  fix(indf_o2)  = 0.2e+6_dp * fakt
  fix(indf_co2) = 400.0_dp * fakt
#ENDINLINE
```

Note that water vapor is considered as fixed as the water vapor variable `q` is computed somewhere else in PALM. In the absence of a prognostic water vapor variable a constant value of q=0.01 kg/kg is assumed.

The unit of `fix` is molecules cm<sup>-3</sup>. `qvap` is already converted to molecules cm<sup>-3</sup>, for O<sub>2</sub> and CO<sub>2</sub> the multiplication by `fakt` converts ppm to molecules cm<sup>-3</sup> .

