# Pollen Emission Model (EMPOL)

## Overview

The pollen emission model implemented in PALM-4U is adopted from EMPOL 1.0 ([Zink et al., 2013](https://doi.org/10.5194/gmd-6-1961-2013)). The EMPOL model for PALM-4U includes four pollen species as in the the original model EMPOL, namely Birch tree, Alder tree, Grasses and Ambrosia. The parameterization of meteorological influence on Ambrosia is still undergoing testing, therefore, Ambrosia has not been activated as yet.  Development of plant and release of the pollen are the processes  specific to plant species and varies significantly, therefore EMPOL parameterization for PALM-4U has been kept flexible with regard to biological and physical processes so that it can be adopted to different plant species.

The emission parameterization in the EMPOL model takes into account the time and amount of emissions of pollen of a given plant species. The former is based on the seasonal model by ([Helbig et al., 2004](https://doi.org/10.1023/B:AERO.0000022984.51588.30)), and later is based on the  ([Zink et al., 2013](https://doi.org/10.5194/gmd-6-1961-2013)). Currently implementation of EMPOL into PALM describes only the  emission  parameterization of pollen grains which is a combination of physical and  biological processes that are characteristic of each plant species. For seasonal description, the user is given the control to define emission time of a given pollen species. The complete seasonal model will be added to PALM-4U as a separate utility code in a future version of the pollen model

The emission process is separated into two steps. The first step is the release of the pollen from the flowers(anthesis) into a pollen reservoir (pollen presentation). This is driven by temperature and relative humidity. The second step is the entrainment of the pollen from the reservoir into the atmosphere provided the meteorological conditions are favorable. This is driven by the turbulent kinetic energy. In some cases when local wind turbulence is intense,  pollen grains can also be released from the flowers and entrained directly into the atmosphere. The emission parameterization of EMPOL comprised of three integral parts a) description of the pollen season adopted from seasonal model from ~\cite{helbig2004numerical}, b) pollen emission and , c) pollen dispersion which includes transport, deposition, sedimentation and washout by precipitation of pollen grains. The pollen parameterization adopted in PALM-4U comprised of the later two processes. The user is given the choice to describe the state of the pollen season state, see section  <a href="#T01">NAMELIST INPUT PARAMETERS</a> for further details.

The pollen parameterization can be described as the maximal daily amount of pollen  $Q_{(pol\_day)}$ that could be released at the peak of the pollen season on one square meter if the conditions for pollen release and entrainment were perfect. The maximum release of a particular pollen species can be tuned by calculating the overall bias between measurements and model. The emitted pollen is then the adjusted amount of $Q_{(pol\_day)}$ by various biological and physical factors. All pollen parameterization factors takes values between 0 and 1 and describe resistence to the pollen release ([Zink et al., 2013](https://doi.org/10.5194/gmd-6-1961-2013)). The amount of pollen that can be released per time step  is calculated considering the following requirements/assumptions:

a) Under optimal conditions for pollen release, the flowers can run out of pollen grains before the end of the day.
b) The daily cycle of pollen release is not prescribed in the model. The amount of released pollen does not depend on the time of the day but results from the temperature and relative humidity factors.

The main assumption in the pollen parameterization is that pollen are dispersed instantaneously and homogeneously. Since processes leading to pollen release are slightly different for different plant species, therefore pollen release functions are plant-dependent and need to be adapted for the different species, whereas processes related to entrainment are driven by meteorological conditions namely moisture and turbulence. The pollen emission flux equation $ F_{(E,pol)}$ of the two (release and entrainment) processes is given by:
$$
 F_{(E,pol)} = \left[\frac{R_{(pol)} } {\Delta t} \times f_{(E,TKE)} \times f_{(E,RH)}
\right]
$$

where the first term to the right $\frac{R_{(pol)} } {\Delta t}$ describes the total net release of pollen per time-step from flower to the reservoir, the second term $f_{(E,TKE)}$ describes the entrainment by turbulent winds, where as the third term $f_{(E,RH)}$ describes the changes in the amount of pollen entrainment into the atmosphere due to moisture.


## <b id="T01"> NAMELIST  (<code>_p3d</code>)   INPUT PARAMETERS </b>

Following is the list of input parameters that can be used to steer EMPOL with PALM-4U.  All namelist option must be listed in the (<code> &chemistry_parameters</code>).


<table style="width:98%;">
<colgroup>
<col style="width: 24%" />
<col style="width: 3%" />
<col style="width: 5%" />
<col style="width: 66%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Parameter name</th>
<th style="text-align: left;">Fortran type</th>
<th style="text-align: left;">Default value</th>
<th style="text-align: center;">Explanation</th>
</tr>
</thead>
<tbody>

<tr class="odd">
<td style="text-align: left;">emis_pollen</td>
<td style="text-align: center;">L</td>
<td style="text-align: center;">.F.</td>
<td style="text-align: left;"> Parameter to switch pollen emissions ON (.TRUE.) or OFF (.FALSE.).</td>
</tr>
<tr class="even">
<td style="text-align: left;">epol_ignore_precip</td>
<td style="text-align: center;">L</td>
<td style="text-align: center;">.F.</td>
<td style="text-align: left;">Parameter to switch precipitation in OFF (.TRUE.) or include  (.FALSE.). in the pollen release processes. Currently precipitation is turned off in the pollen model code.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epol_ignore_solar</td>
<td style="text-align: center;">L</td>
<td style="text-align: center;">.F.</td>
<td style="text-align: left;">Switch to ignore (.TRUE.) or include solar radiation (.FALSE.) in the pollen release prrocesses. Currently soloar activity is turned off in the pollen model code.</td>
</tr>
<tr class="even">
<td style="text-align: left;">epol_model</td>
<td style="text-align: center;">C</td>
<td style="text-align: center;">'zink'</td>
<td style="text-align: left;"> This parameter offers the choice of pollen emissin parameterizations. Currently only 'zink' <a href="https://gmd.copernicus.org/articles/6/1961/2013/">(Zink et al., 2013) </a>, is implemented and available. In the future releases of the pollen model, 'hlebig' <a href="https://doi.org/10.1023/B:AERO.0000022984.51588.30">(Helbig et al., 2004) </a>, and 'sofiev', <a href="https://doi.org/10.1007/s00484-012-0532-z">(Sofiev et al., 2013) </a>, would be included.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epol_pool_reset_hour</td>
<td style="text-align: center;">I</td>
<td style="text-align: center;">0</td>
<td style="text-align: left;">The time in UTC when the pollen emission pool is reset and replenished. In most of the cases, it is 00:00 hours local time.</td>
</tr>
<tr class="even">
<td style="text-align: left;">epol_seasonal_factors</td>
<td style="text-align: center;">R</td>
<td style="text-align: center;">0.9</td>
<td style="text-align: left;">Seasonal description of pollen species. The parameter varies between 0 and 1. This parameter describes start, end and peak of the pollen season. A value of 0 means no pollen emission season, and a value close to 1 indicates peak of pollen emission season.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epol_specs_names</td>
<td style="text-align: center;">C</td>
<td style="text-align: center;">'no_value'</td>
<td style="text-align: left;">Names of the pollen emitting plant species. Currently thee species are available, namely, birch (Betula), Alder, and grasses (Poaceae). Due to reasons for algorithm testing and authentication, parameterization for Ambrosia has not not been activated as yet. Ambrosia will be availabe in the next release of the pollen model.</td>
</tr>
<tr class="even">
<td style="text-align: left;">epol_tke_scheme</td>
<td style="text-align: center;">C</td>
<td style="text-align: center;">'default'</td>
<td style="text-align: left;">Parameterization to estimate TKE. Currently three TKE schemes are available namely, 'default' 'dynamic', and 'adhoc'.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epol_tke_sgs_fraction</td>
<td style="text-align: center;">R</td>
<td style="text-align: center;">0.1</td>
<td style="text-align: left;">The TKE fraction is required to calculate TKE for the 'default' TKE scheme 'epol_tke_scheme'. The value may vary between > 0 and 1. The greater the fraction, the higher the magnitude of the TKE.</td>
</tr>
<tr class="even">
<td style="text-align: left;">epol_tree_specs</td>
<td style="text-align: center;">I</td>
<td style="text-align: center;">999</td>
<td style="text-align: left;">Pollen emitting tree types. The current PALM database has 86
<a href="https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids/palm_csd#Singletrees">tree types </a>. The EMPOL model calculates pollen from Birch (betula) and Alder trees. The user can choose other pollen emitting trees also from this database, however, the user need to provide the pollen attributes for the particular tree pollen.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epol_tuning_factors</td>
<td style="text-align: center;">R</td>
<td style="text-align: center;">1.0</td>
<td style="text-align: left;">This parameter is used to tune the maximum pollen release per day m<sup>2</sup> of a particular pollen species based on the observations. The value must be greater than 0. There is no upper bound for this parameter, however, in most of the cases it will remain less than 10.</td>
</tr>
<tr class="even">
<td style="text-align: left;">epol_update_interval</td>
<td style="text-align: center;">R</td>
<td style="text-align: center;">300.0</td>
<td style="text-align: left;"> Time interval to update the pollen model. In case the user defined update interval is less than model's physical time-step, then model will issue a warning and set the pollen model time-step to PALM model time-step.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">epol_vegetation_spcs</td>
<td style="text-align: center;">I</td>
<td style="text-align: center;">999</td>
<td style="text-align: left;">vegetation type. PALM utilizes ECMWF-IFS classification for <a href="https://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#vegetation_type">vegetation types</a>  The EMPOL model calculates pollen from grasses(poaceae) and Ambrosia(Ragweed). Currently only grasses (Poaceae) vegetation type is available. The user can also choose other pollen emitting vegetation types from this database, however, the user need to provide the pollen attributes for the selected vegetation. type</td>
</tr>
</tbody>
</table>

## EMISSION INPUT
The pollen model does not require a separate emission input file. However, the geostatic (<code>_static</code>) file must contain two fields, namely, <code>vegetation_type</code> (for grasses and 2D vegetation patches) and <code>tree_type</code> (for 3D resolved trees). 

## LIMITATION
1) The present implementation of EMPOL in PALM-4U only considers three pollen species, namely, Birch tree, Alder tree, and Grasses. The forth pollen species i.e Ambrosia will be added later.
2) Seasonal description model is currently not added to the pollen model. A namelist parameter, <code>epol_seasonal_factors</code> is provided to set the emission time of pollen species by the user based on measurements.
3) Altitude effects pollen emission of trees. Given the relatively smaller urban PALM domains, currently, altitude is not considered and set to 1.
4) Grid cell fraction of plants is not considered in this version of EMPOL.
5) The user is responsible to provide vegetation information in the format as described on [static driver](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids/static) page.

## EXAMPLE SETUP
Following is the example chemistry namelist for steering pollen model. 

```
!--------------------------------------------------------------------------------
!-- CHEMISTRY MODEL PARAMETER NAMELIST
!   Documentation: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/chempar
!
!   Chemsistry pollen namelist for 20th April (peak emission day of Birch trees)
!--------------------------------------------------------------------------------
 &chemistry_parameters

    chem_gasphase_on          = .TRUE.,
    chem_mechanism            = "empol1.0",
    photolysis_scheme         = 'simple',
    emissions_anthropogenic   = .FALSE.,
    call_chem_at_all_substeps = .FALSE.,
    bc_cs_b                   = 'neumann',
    bc_cs_t                   = 'neumann',
    emis_pollen               = .TRUE.,
    epol_ignore_precip        = .TRUE.,
    epol_ignore_solar         = .TRUE.
    epol_update_interval      =  60.0,
    epol_pool_reset_hour      =  0,
    epol_model                = 'zink',
    epol_tke_scheme           = 'default',             !< 'dynamic', 'adhoc',
    epol_tke_sgs_fraction     = 0.10,
    epol_specs_names          = 'POL_ALNU', 'POL_POAC', 'POL_BETU',
    epol_seasonal_factors     =   0.0,         0.0,         0.9,    ! Peak pollen release period:
                                                                    ! Birch   = 0.9 : 15 Apr - 25 Apr.
                                                                    ! Grasses = 0.9 : 15-May - 15 Jul.
                                                                    ! Alder   = 0.9 : 20 Feb - 15 Mar.
                                                                                                 
    epol_tuning_factors       =   999,        0.01,          10,
    epol_tree_specs           =   999,         999,           7,
    epol_vegetation_specs     =   999,           3,         999,

    icntrl(3)                 = 1,    !solver ros2
    icntrl(4)                 = 500,  !max number of chem-substeps
    rcntrl(3)                 = 0.1,  ! Hstart, starting value for the integration step size

 / ! end of chemistry_parameters namelist
```

## CHEMISTRY MECHANISM
A default pollen mechanism (<code>empol1.0</code>) has been added to  readily available chemistry mechanisms, The passive chemistry mechanism contains the same four pollen species to steer the pollen model. User may add/create pollen mechanisms with the reduced/increased number of same/different pollen species. User may combine pollen speceis with other gas-phase and/or aerosol chemistry mechanisms. Detailed instructions as how to create a chemistry mechanism are available on [chemistry mechanisms](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/chemmech) page.


## REFERENCES
Helbig, N., Vogel, B., Vogel, H., & Fiedler, F. (2004). Numerical modelling of pollen dispersion on the regional scale. Aerobiologia, 20(1), 3–19. https://doi.org/10.1023/B:AERO.0000022984.51588.30

Sofiev, M., Siljamo, P., Ranta, H., Linkosalo, T., Jaeger, S., Rasmussen, A., Rantio-Lehtimaki, A., Severova, E., & Kukkonen, J. (2013). A numerical model of birch pollen emission and dispersion in the atmosphere. Description of the emission module. International Journal of Biometeorology, 57(1), 45–58. https://doi.org/10.1007/s00484-012-0532-z

Zink, K., Pauling, A., Rotach, M. W., Vogel, H., Kaufmann, P., & Clot, B. (2013). EMPOL 1.0: A new parameterization of pollen emission in numerical weather prediction models. Geoscientific Model Development, 6(6), 1961–1975. https://doi.org/10.5194/gmd-6-1961-2013
