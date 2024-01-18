---
title: Overview
---
# FASTv8 Coupler

---

!!! warning
    This site is  Work in Progress.

## Coupling Concept

FAST is an open source software package that enables the precise simulation of the various dynamic processes and physical subsystems of a wind turbine. In addition to the ElastoDyn and ServoDyn modules, which handle the simulation of loads and deformations on the rotor and gearbox on the one hand and control and power generation on the other, the aerodynamics of the wind turbine are simulated with the AeroDyn module using e.g. the blade element momentum (BEM) method. BEM is an established method for the calculation of local forces on propellers or wind turbine blades. In most cases, pre-calculated wind fields are used and simplified atmospheric conditions are assumed. The module InflowWind feeds AeroDyn with these flow data.

PALM already contains an internal Actuator Disc Model with Rotation (ADM-R). This is a simple parameterization of a wind turbine without the resolution of individual rotor blades and is specifically tailored to a particular turbine. The coupling of PALM and FAST described here allows the modeling of both the wind turbine for load calculation and the atmospheric conditions with the respective specialized tool. In order to consider also the influence of the wake in larger wind farms, a bi-directional coupling is implemented. A coupling via an Actuator Line Model (ALM) would be the most realistic approach. Here, each blade is modeled by a line divided into blade element points. There is no additional averaging or smearing. However, with bi-directional coupling, PALM would be constrained to the comparatively very small time step of FAST. Also spatially, the rotor would have to be sufficiently high resolution to actually accommodate the information density of ALM in the flow. Especially for larger simulations of e.g. whole wind farms, a high spatial resolution of the LES is often limited by the available computer resources.

The implemented coupling approach therefore uses to the concept of an Actuator Segment Model (ASM) and is based on the work of [Krüger et al. (2022)](https://doi.org/10.5194/wes-7-323-2022). In the 3D flow, the rotor acts via segments that correspond to the circular section swept per time step. The larger the PALM time step compared to the FAST time step, the larger the area of the circular segment that the rotor sweeps during a PALM time step. The rotor segment thus covered serves as a projection surface of the acting forces. This is a computationally efficient middle ground where both PALM and FAST can use the respective internal time step.

Multiple different turbines can be places into the PALM setup simultanously. Each individual turbine runs its own FAST instance and PALM communicates with all these instances. The interface used is generic with no direct software dependency. Sockets and TCP/IP were chosen as a simple coupling with easily exchangeable receivers. All components are addressed by their network address (IP and port). The communication is done by data packed in individual custom messages over the specified TCP sockets.


## FASTv8 Installation

In order use the PALM to FASTv8 coupling, you need to build a modified version of FASTv8 first. You can get the modified FASTv8 version [here](https://gitlab.palm-model.org/fast/fastv8). Please follow the instructions in Compiling/CompilingInstructions_FASTv8.pdf in order to build FASTv8 (see page 11 and 12 if you are working on Linux).

For a scientific and technical introduction to the Coupler, please see [Krüger et al. (2022)](https://doi.org/10.5194/wes-7-323-2022). In Short, the two Models PALM and FASTv8 both communicate using a custom TCP/IP interface. They requrie a matching setup as well as information about host IP and port number on which the respective other model is listening.

PALM also requires some small modifications during configuration and installation in order to use the PALM tp FASTv8 coupler. The cpp macro `__fastv8` needs to be activated. This can be done by adding `-D__fastv8` to the variable `cpp_options` in your PALM config file. Additionally, the variables `c_compiler_name` and `c_compiler_options` must be set according to your compiler in your PALM config file.

## Usage

The modified version of FASTv8 contains a example Setup called NRELOffshrBsline5MW. The following usage example is based on this Setup.

Each FAST instance receives an individual ID and is also told at which network address it should wait for PALM to connect. The TCP/IP connection configuration is done inside the file `palm_sr.srv` and contains the following lines:

```
TURBINE_ID 1
HOST_FAST localhost
PORT_FAST 20001
```

PALM receives a list with the network addresses of all FAST instances to be expected. Here the TCP/IP connection configuration is done inside the _p3d namelist `fastv8_coupler_parameters` (see below).

First all FAST instances are started (they initally run idle). In order to run this example, please enter the directory `NRELOffshrBsline5MW` inside the directory of the modified FASTv8 version and execute the following command:
```bash
../bin/FAST_glin64 NRELOffshrBsline5MW_Onshore.fst -c
```

Then PALM is started using for example the following minimal example setup:

```fortran
!-------------------------------------------------------------------------------
!-- INITIALIZATION PARAMETER NAMELIST
!   Documentation: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/inipar
!-------------------------------------------------------------------------------
&initialization_parameters
!
!-- grid parameters
!-------------------------------------------------------------------------------
    nx                         = 47, ! Number of gridboxes in x-direction (nx+1)
    ny                         = 95, ! Number of gridboxes in y-direction (ny+1)
    nz                         = 48, ! Number of gridboxes in z-direction (nz)

    dx                         = 10.0, ! Size of single gridbox in x-direction
    dy                         = 10.0, ! Size of single gridbox in y-direction
    dz                         = 10.0, ! Size of single gridbox in z-direction
!
!-- initialization
!-------------------------------------------------------------------------------
    initializing_actions       = 'set_constant_profiles', ! initial conditions
    reference_state            = 'horizontal_average',

    damp_level_1d              = 1000.0,
    dt_run_control_1d          = 86400.0,
    dt_pr_1d                   = 864000.0,

    ug_surface                 = 8.0, ! u-comp of geostrophic wind at surface
    vg_surface                 = 0.0, ! v-comp of geostrophic wind at surface

    omega                      = 0.0,

    pt_surface                 = 283.15, ! initial surface potential temp

    pt_vertical_gradient       =    0.0,
    pt_vertical_gradient_level =    0.0,


!
!-- boundary conditions
!-------------------------------------------------------------------------------
    constant_flux_layer        = .FALSE.,
    bc_pt_b                    = 'neumann', ! required with surface_heatflux
    bc_uv_b                    = 'neumann',
    bc_lr                      = 'dirichlet/radiation',
!
!-- numerics
!-------------------------------------------------------------------------------
    psolver                    = 'multigrid',

/ ! end of initialization parameter namelist

!-------------------------------------------------------------------------------
!-- RUNTIME PARAMETER NAMELIST
!   Documentation: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/d3par
!-------------------------------------------------------------------------------
&runtime_parameters
!
!-- run steering
!-------------------------------------------------------------------------------
!   dt                         = 0.5
    end_time                   = 60.1, ! simulation time of the 3D model

    create_disturbances        = .TRUE., ! randomly perturbate horiz. velocity

!
!-- data output
!-------------------------------------------------------------------------------
    netcdf_data_format         = 5, ! use NetCDF4 with parallel IO
    interpolate_to_grid_center = .TRUE.,

    dt_run_control             = 0.0,    ! output interval for run control
    dt_do3d                    = 1.12,


    data_output                = 'u', 'v', 'w', 'ti', 'thrx', 'tory', 'torz',

/ ! end of runtime parameter namelist

!-------------------------------------------------------------------------------
!-- FASTv8 COUPLER PARAMETER NAMELIST
!-------------------------------------------------------------------------------
&fastv8_coupler_parameters
!
!-- turbine steering
!-------------------------------------------------------------------------------
    time_turbine_on            = 0.0,

    nturbines                  =   1,

    fast_n_blades_max          =   3,
    fast_n_blade_elem_max      =  62,
    fast_host_addr_default     = '127.0.0.1',

!    fast_host_addr(1)          = '127.0.0.1',
    fast_host_port(1)          = '20001',
    palm_tower_ref_pos_x(1)    =  350.0,
    palm_tower_ref_pos_y(1)    =  150.0,
    palm_tower_ref_pos_z(1)    =    0.0,
    turb_C_d_tow(1)            =    1.2,
    dtow(1)                    =    3.87,
    htow(1)                    =    87.6,

    fast_host_port(2)          = '20002',
    palm_tower_ref_pos_x(2)    =  350.0,
    palm_tower_ref_pos_y(2)    =  350.0,
    palm_tower_ref_pos_z(2)    =    0.0,

    fast_host_port(3)          = '20003',
    palm_tower_ref_pos_x(3)    =  350.0,
    palm_tower_ref_pos_y(3)    =  550.0,
    palm_tower_ref_pos_z(3)    =    0.0,

    fast_host_port(4)          = '20004',
    palm_tower_ref_pos_x(4)    =  350.0,
    palm_tower_ref_pos_y(4)    =  750.0,
    palm_tower_ref_pos_z(4)    =    0.0,

!    reg_fac                    = 2.0,
!    fbox_fac                   = 1.5,

/ ! end of fastv8 coupler parameter namelist

```

After all connections are establish, the coupled simulation starts executing. In a multi-turbine setup, each turbine can be operated with individual location, rotor radius, hub heights, power characteristics and turbine controller (The full set of FAST configuration files can be used). PALM gets the detailed properties of each turbine directly from the corresponding FAST instance.

## Notes, shortcommings and open issues

1. The ASM approach can show its strengths when the PALM time step is a higher multiple of the FAST time step. If the difference is too large, so that the rotor sweeps the entire rotor plane within one PALM time step, then the use of ADM is more appropriate. An application with very high spatial resolution of the turbines in the PALM grid is possible, but the ASM approach then behaves more like an ALM without being evaluated for it. The averaging and smearing used is then a likely source of error.

2. The current implementation is based on FASTv3, but the current successor of FASTv8 is a rewrite called [OpenFAST](https://github.com/OpenFAST/openfast). We recommend porting this coupler to the newest OpenFAST version before any continued Development.

3. Further optimization and parallelization is also urgently required. As far as possible, the focus should be on asynchronous execution of the two components. Currently the communication is blocking.
