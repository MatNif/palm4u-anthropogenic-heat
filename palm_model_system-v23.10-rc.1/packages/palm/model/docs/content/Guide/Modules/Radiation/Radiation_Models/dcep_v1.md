# Urban canopy parametrization scheme for coarse-grid simulation 

When simulating large areas of a city using PALM in a coarse-grid resolution (e.g. ~ 100m), buildings and obstacles can not be resolved explicitly. However, the dynamic and thermodynamic effects of buildings and obstacles are implicitly considered in these non-building resolving simulations to represent the interaction of atmospheric flow with urban areas. This consideration is realized using the urban Double-Canyon Effects Parametrization scheme (DCEP) [Schubert et al. 2012](https://doi.org/10.1007/s10546-012-9728-3).

## DCEP scheme

DCEP is a multi-layer urban canopy parametrization scheme [Schubert et al. 2012](https://doi.org/10.1007/s10546-012-9728-3) based on the Building Effect Parametrization (BEB) [Martilli et al. 2002](https://doi.org/10.1023/A:1016099921195). DCEP scheme calculates the incoming and outgoing longwave and shortwave radiation for the surfaces of an urban street canyon, i.e. the roof, wall and ground surfaces. The building morphology of the urban area is used to characterize the urban street canyon using its street and building width as well as its canyon length, and the height distribution of buildings. 

The DCEP scheme is described in [Schubert 2013](https://refubium.fu-berlin.de/handle/fub188/11568) and [Schubert et al. 2012](https://doi.org/10.1007/s10546-012-9728-3) and was originally coupled with the mesoscale model COSMO/CLM [Schubert and Grossman 2013](https://doi.org/10.1002/qj.2311).

## Coupling DCEP scheme to PALM model system

DCEP scheme is coupled online with the PALM model system. The DCEP scheme is located at a separate module, i.e. `dcep_mod.f90`, and all the related subroutines are named so that they have the prefix `dcep_`. In the initialization phase, DCEP calculates the view factors, turbulent length scales, and the radiation matrices. In the time integration phase (time stepping), DCEP receives the the incoming longwave and shortwave radiation fluxes from the radiation scheme used in PALM, e.g. RRTMG, and calculates the urban radiation fluxes, based on the view factors. These radiation fluxes are used to calculate the heat storage and the fluxes from buildings, based on the interpolation of the dynamic quantities. The fluxes from buildings are aggregated to the coarse PALM grid to calculate the surface fluxes, turbulent fluxes, and the dynamics on the PALM grid.

## Using DCEP/PALM coupled system

### Installation

Compiling the DCEP/PALM model system does not require adding extra special C pre-processor directives (`%cpp_options`) to the palm configuration file. However, if the parallel Linear Algebra Package (PLAPACK) is to be used for matrices multiplication, the string `-D__DCEPLAPACK` has to be added to `%cpp_options` field.

### Input Data

DCEP/PALM model system requires the information about urban area such as horizontal building size, street canyon width and a building density as a function of height. These data should be provided via a NetCDF file and the file name should be `<run_identifier>_dcep`. The variables are listed in the following table:

| Variable name | Description | Size | Unit |
| --- | --- | --- | --- |
| BUILD_PROP | probability to have a building at the height | (N<sub>uc</sub>, N<sub>fdir</sub>, N<sub>uh</sub>, ny, nx) | - |
| BUILD_W | building width in grid cell | (N<sub>uc</sub>, N<sub>fdir</sub>, ny, nx) | m |
| FR_BUILD | fraction of building surface in grid cell | (N<sub>uc</sub>, ny, nx) | - |
| FR_STREETD | street fraction | (N<sub>uc</sub>, N<sub>fdir</sub>, ny, nx) | - |
| FR_URBAN | fraction of urban surfaces in grid cell | (N<sub>uc</sub>, ny, nx) | - |
| FR_URBANCL | urban classes fraction | (N<sub>uc</sub>, ny, nx) | - |
| STREET_W | street width in grid cell | (N<sub>uc</sub>, N<sub>fdir</sub>, N<sub>uh</sub>, ny, nx) | m |

The input/output configuration file, `.palm.iofiles`, should be adjusted to consider the input file `<run_identifier>_dcep` needed to run DCEP/PALM model as follows:

``` bash
PIDS_DCEP inopt:tr d3#:d3r $base_data/$run_identifier/INPUT _dcep
```
### Changes in namelists

#### DCEP parameters namelist

To activiate DCEP, an extra namelist should be added to the control file `<run_identifier>_p3d`, i.e. `&dcep_parameters`. This namelist contains all the settings and options needed to steer the DCEP scheme. The following table contains all the possible variables.

| Variable name | FORTRAN Type | Default | Description | unit |
| --- | --- | --- | --- | --- |
| albedo_ground | R | 0.2 | albedo of the ground | - |
| albedo_roof | R | 0.2 | albedo of the roof | - |
| albedo_wall | R | 0.2 | albedo of the wall | - |
| dcep_average_radiation | L | .F. | flag to activiate average_radiation from DCEP module | - |
| dzlayer_ground | R | -1.0 | thickness of ground's layers | m |
| dzlayer_roof | R | -1.0 | thickness of roof's layers | m |
| dzlayer_wall | R | -1.0 | thickness of wall's layers | m |
| emissivity_ground | R | 0.95 | emissivity of ground | - |
| emissivity_roof | R | 0.9 | emissivity of roof | - |
| emissivity_wall | R | 0.9 | emissivity of wall | - |
| iurb_cdrag | I | 1 |  type of drag coefficient for walls | - |
| iurb_ls | I | 2 | type of urban length scale parametrization | - |
| ke_ground | I | 10 | number of grid levels in the ground | - |
| ke_roof | I | 10 | number of grid levels in roofs | - |
| ke_uhl | I | 13 | number of urban height levels | - |
| ke_wall | I | 10 | number of grid levels in walls | - |
| limpufl | L | .T. | treat urban tendencied implicitly flag | - |
| lrroofs | L | .F. | flag to use roof radiation for budget | - |
| ltintfix | L | .F. | fixed innermost temperature in urban surfaces | - |
| lurbradcor | L | .T. | use correction factor for radiation reduction factor for radiation from canyon sides | - |
| lurbvel | L | .T. | urban modification of wind velocity | - |
| n_uclass | I | 1 | number of urban classes in DCEP model | - |
| n_udir | I | 4 | maximum number of street directions | - |
| rlength_ground | R | 0.01 | ground's roughness length | m |
| rlength_roof | R | 0.01 | roof's roughness length | m |
| thermdiff_ground | R | 0.67E-6 | thermal diffusivity of ground's layers | - |
| thermdiff_roof | R | 0.29E-6 | thermal diffusivity of roof's layers | - |
| thermdiff_wall | R | 0.67E-6 | thermal diffusivity of wall's layers | - |
| tinterior_ground | R | 292.0 | ground's initial and interior temperature | K |
| tinterior_roof | R | 292.0 | roof's initial and interior temperature | K |
| tinterior_wall | R | 292.0 | walls's initial and interior temperature | K |
| z_uhl | R | -99999.9 | heigth of ith face of level in urban class | m |

#### Runtime parameters namelist

The output parameters related to DCEP can be added to the output list data_output `&runtime_parameters` in the control file `<run_identifier>_p3d`. The name of all DCEP related outputs start with `dcep_` to be identified. A list of the currently available output parameters are listed in the table below.

| Variable name | Description | unit |
| --- | --- | --- | 
| dcep_nz_meso* | number of mesoscale height levels in urban layer | - |
| dcep_fr_wall | Probability that a building has an height greater or equal to z, corresponds to roof heigths | - |
| dcep_fr_roof | probability that a building has an height equal to z | - |
| dcep_tt_urb | tendency of temperature change due to urban area | K s -1 |
| dcep_shfl_roof | sensible heat flux from roofs | W m -2 |
| dcep_sw_roof | incoming shortwave radiation on roofs | W m -2 |
| dcep_t_g_urb* | effective urban ground temperature | K |
| dcep_shfl_urb | total sensible heat flux from urban surfaces | W m -2 |
| dcep_albedo_urb* | effective urban albedo | - |
| dcep_emiss_urb* | urban emissivity | - |
| dcep_t_grad_urb* | effective urban radiation temperature | K |
| dcep_rl_roof | incoming longwave radiation on roofs | W m -2 |
| dcep_rl_wallw | incoming longwave radiation on west wall | W m -2 |
| dcep_rl_walle | incoming longwave radiation on east wall | W m -2 |
| dcep_strfl_urb* | total ground flux from urban surfaces | W m -2 |
| dcep_t_ground | temperature in each layer of the ground | K |
| dcep_t_roof1 | temperature at first layer of the roof | K |
| dcep_t_walle | temperature at first layer of the east wall | K |
| dcep_t_wallw | temperature at first layer of the west wall | K |
