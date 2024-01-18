# ISORROPIA model for SIA

## Overview

For the modelling of production of secondary inorganic aerosols (SIA), the user can choose to use subroutines directly from the third-party aerosol thermodynamic equilibrium models ISORROPIA (Nenes, 1998) and ISORROPIA II (Fountoukis and Nenes, 2007). For convenience, both ISORROPIA and ISORROPIA II will be referred to collectively as ISORROPIA, except in cases where distinction is necessary.  Further information pertaining to the application of ISORROPIA can be found, for instance, in the ISORROPIA Reference Manual (Fountoukis et al, 2009).

The use of ISORROPIA in PALM can be activated from the namelist (<code>_p3d</code>  file), where the ISORROPIA solver behavior can be further defined.  Source codes of ISORROPIA and ISORROPIA II can be requested directly from the respective authors.  While these options are designed to reflect the functionality of different aspects of the ISORROPIA solvers, some terminologies specific to ISORROPIA will be used for those who are already familiar with ISORROPIA. The present documentation will provide build instructions for ISORROPIA in PALM and to provide namelist options for activating and controlling the behavior of the ISORROPIA solver.

## Build Instructions

ISORROPIA must be compiled as a standalone shared library (<code>.so</code>) with the same   toolchain used in building PALM.  The make and version of the toolchain can be verified by inspecting the corresponding PALM build configuration file (<code>.palm.config.*</code>) under the entry <code>%compiler_name%</code> or <code>%compiler_name_ser%</code>.   The code must be compiled as position-independent code into a shared object file, which can be accomplished respectively with the build switches <code>-fPIC</code> and <code>-shared</code> using the GNU or Intel compiler toolchain, for example.

Once the <code>.so</code> file has been successfully built, the PALM configuration must be modified so that the ISORROPIA library can be integrated, in which corresponding changes must be made in the entries for <code>%compiler_options</code> and <code>%linker_options</code> in the <code>.palm.config</code> file. The following table outlines the necessary changes.

| Entry in <code>.palm.config</code> | Additional Option |
|---|---|
| <code>%compiler_options</code> | <code>-D__ISORROPIA1</code> (ISORROPIA) *or* <br><code>-D__ISORROPIA2</code> (ISORROPIA II) |
|<code>%linker_options</code>          | Path to the ISORROPIA <code>.so</code> file |

Once these changes have been made, <code>palmbuild</code> can be executed to include the ISORROPIA model.  Please note that for each build either ISORROPIA or ISORROPIA II, but not both, can be integrated into PALM. 
 
## Namelist (<code>_p3d</code>) Options
 
To activate the use ISORROPIA library subroutines in PALM, the namelist option <code>chem_isorropia</code> is to be set to <code>.TRUE.</code>. Optionally, the frequency at which the ISORROPIA solver is called can be set using the option <code>chem_isorropia_update_interval</code>, where a numerical value given indicates the simulation time in seconds between successive calls to the ISORROPIA solver. The default is set to 300 seconds. 

In addition, while the ISORROPIA solver will operating using default parameters values, they can be customized from the namelist.  These options are presented in the table below, along with the corresponding variable names in ISORROPIA as a reference for users who are familiar with ISORROPIA.  Please note that all namelist items pertaining to ISORROPIA will begin with the prefix <code>chem_isorropia</code>. Detailed explanation of the different options can be found in the ISORROPIA Reference Manual (Fountoukis et al, 2009).  Any invalid parameter input in the namelist will be replaced with the corresponding default value.

| Namelist (<code>_p3d</code>) Option (<code>chem_isorropia_</code>) | ISORROPIA Variable | Function | Possible Values (Default in Boldface) | 
|---|---|---|---|
| <code>_problem_type</code> | <code>CNTRL(1)</code> | Solver mode | **0.0**, 1.0 |
| <code>_aerosol_state</code> | <code>CNTRL(2)</code> | Aerosol state | **0.0**, 1.0 |
| <code>_max_iteration</code> | <code>MAXITI</code> | Maximum solver iterations | Any positive integer (**100**) |
| <code>_solver_tolerance</code> | <code>EPSI</code> | Solver convergence criterion | Any positive floating point (**1.0E-6**) |
| <code>_activity_coefficient_method</code> | <code>IACALCI</code> | Activity coefficients algorithm | 0, 1 |
| <code>_max_activity_sweep</code> | <code>NSWEEPI</code> | Maximum activity coefficient sweeps  | Any positive integer (**4**) |
| <code>_activity_tolerance</code> | <code>EPSACTI</code> | Activity coefficient convergence criterion | Any positive floating point (**0.05**) |
| <code>_root_subdivisions</code> | <code>NDIVI</code> | Number of subdivisions for root tracking | Any positive integer (**5**) |
| <code>_mdr_weight_method</code> | <code>WFTYPI</code> | Type of weighting algorithm for mutual deliquescence regions (MDR) | **0**, 1, 2 |
| <code>_mass_conservation_mode</code> | <code>NADJI</code> | Mass conservation enforcement (ISORROPIA II only) | 0, **1** |

## Aerosol Species in Chemical Mechanism

PALM will collect output aqueous and solid aerosol species which are defined in the active chemical mechanism.  The following table lists the aerosol species that are available form ISORROPIA.  Please note that not all species listed below needed be included in the PALM chemical mechanism; only those that are required.  In addition, computation of aerosol species involving Magnesium (Mg), potassium (K), and calcium (Ca) are only available in ISORROPIA II.

| Aqueous Aerosol | PALM Mechanism Name |  |  Solid Aerosol | PALM Mechanism Name |
|---|---|---|---|---|
| H<sup>+</sup> | H_aq | | NaNO<sub>3</sub> | NANO3_s |
| Na<sup>+</sup> | NA_aq | | NH<sub>4</sub>NO<sub>3</sub> | NH4NO3_s |
| NH<sub>4</sub><sup>+</sup> | NH4_aq | | NaCl | NACL_s |
| Cl<sup>-</sup> | CL_aq | | NH<sub>4</sub>Cl | NH4CL_s |
| SO<sub>4</sub><sup>2-</sup> |SO4_aq | | Na<sub>2</sub>SO<sub>4</sub> | NA2SO4_s |
| HSO<sub>4</sub><sup>-</sup> | HSO4_aq | | (NH<sub>4</sub>)<sub>2</sub>SO<sub>4</sub> | NH42SO4_s |
| NO<sub>3</sub><sup>-</sup> | NO3_aq | | NaHSO<sub>4</sub> | NAHSO4_s |
| H<sub>2</sub>O | H2O_aq | | NH<sub>4</sub>HSO<sub>4</sub> | NH4HSO4_s |
| NH<sub>3</sub> | NH3_aq | | (NH<sub>4</sub>)<sub>4</sub>H(SO<sub>4</sub>)<sub>2</sub> | NH44HSO42_s | 
| HCl | HCL_aq | | CaSO<sub>4</sub> | CASO4_s |
| HNO<sub>3</sub> | HNO3_aq | | Ca(NO<sub>3</sub>)<sub>2</sub> | CANO32_s |
| OH<sup>-</sup> | OH_aq | | CaCl<sub>2</sub> | CACL2_s |
| Ca<sup>2+</sup> | CA_aq | | K<sub>2</sub>SO<sub>4</sub> | K2SO4_s |
| K<sup>+</sup> | K_aq | | KHSO<sub>4</sub> | KHSO4_s |
| Mg<sup>2+</sup> | MG_aq | | KNO<sub>3</sub> | KNO3_s |
| | | | KCl | KCL_s |
| | | | MgSO<sub>4</sub> | MGSO4_s |
| | | | Mg(NO<sub>3</sub>)<sub>2</sub> | MGNO32_s |
| | | | MgCl<sub>2</sub> | MGCL2_s |

## References

- Fountoukis, C., and Nenes, A.: ISORROPIA II: a computationally efficient thermodynamic equilibrium model for K+–Ca2+–Mg2+–NH+ 4 –Na+–SO2− 4 –NO− 3 –Cl−–H2O aerosols, Atmospheric Chemistry and Physics, 7, 4639–4659, 2007.
- Fountoukis, C., Nenes, A., Pandis, S., and Pilinis, C.: ISORROPIA v 2.1 Reerence Manual, 2009.
- Nenes, A., Pandis, S. N., and Pilinis, C.: ISORROPIA: A new thermodynamic equilibrium model for multiphase multicomponent inorganic aerosols, Aquatic Geochemistry, 4, 123–152, 1998.
