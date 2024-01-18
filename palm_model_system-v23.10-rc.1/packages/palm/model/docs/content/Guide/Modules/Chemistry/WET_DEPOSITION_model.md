# Wet deposition parameterization

## Overview

Wet deposition in PALM is modelled using the parameterization of Berge and Jakobsen (1998). Scavenging rates for specific soluable gas and particulate species are calculated at the in- and below-cloud levels based on the field precipitation rate, through which concentrations of said species can be updated along each vertical column in the computational domain.

## Theoretical foundations

The convention employed in Simpson et al (2012) is used here to describe the parametrization of scavenging rates, expressed in $\text{s}^{-1}$, which consists of the following three calculations:

1. In-cloud scavenging rate ($S_\text{in}$),
2. Under-cloud scavenging rate for gas phase species ($S_\text{sub}$), and
3. Under-cloud scavenging rate for particulate species ($S_\text{sub}^\text{aer}$).

The in-cloud scavenging rate applies to all relevant gas and particulate species, and is calculated as follows:

$\displaystyle{S_\text{in} = -\left(\frac{W_\text{in}P}{h_S\rho_w}\right)\chi}{\quad\quad}$**(1)**,

where<br>
$W_\text{in}$ [ ] is the in-cloud scavenging ratio,<br>
$P$ [$\text{kg}$ $\text{m}^{-2}$ $\text{s}^{-1}$] is the precipitation rate, <br>
$h_S$  [$\text{m}$] is the characteristic scavenging depth,<br>
$\rho_w$ [$\text{kg}$ $\text{m}^{-3}$] is density of liquid water, and<br>
$\chi$ [ ] is the mixing ratio of a soluble component.

Similarly, the below-cloud scavenging rate for gas-phase species can be modelled using the relation below:

$\displaystyle{S_\text{sub} = -\left(\frac{W_\text{sub}P}{h_S\rho_w}\right)\chi}{\quad\quad}$**(2)**,

where $W_\text{sub}$ [ ] is the in-cloud scavenging ratio.

On the other hand, the scavenging of particle species through entrainment of rain droplets is determined using the following equation:

$\displaystyle{S_\text{sub}^\text{aer} = -\left(\frac{\overline{E}AP}{V_\text{dr}}\right)\chi}{\quad\quad}$**(3)**,
 
where<br>
$\overline{E}$ [ ] is the raindrop particle collection efficiency,<br>
$A (= 5.2 \text{m}^3$ $\text{kg}^{-1}$ $\text{s}^{-1})$ is an empirical constant based on a Marshall-Palmer size distribution (Laakso et al, 2003), and<br>
$V_\text{dr}w$ [$\text{m}$ $\text{s}^{-1}$] is the rain drop terminal velocity.

Typically, the cloud layer is identified by using the mass fraction of liquid water, when it crosses a threshold of 10$^\text{5}$ on a vertical column (van Zanten et al, 2011). Meanwhile, the characteristic scavenging depth ($h_S$) is set to a constant of 1000 m (Berge and Jakobsen, 1998), and a value of 5 $\text{m}$ $\text{s}^{-1}$ has been adapted for the rain drop terminal velocity (Simpson et al, 2012).  In addition, the scavenging ratios ($W_\text{in}$ and $W_\text{sub}$), as well as, the particle collection efficiency ($\overline{E}$) are dependent on the particular soluble gas-phase and particle species.  These are presented in the table below, adapted from Simpson et al (2012). 

| Species                | $W_\text{in} \times 10^{-6}$ | $W_\text{sub} \times 10^{-6}$ | $\overline{E}$ |
|------------------------|------------------------------|-------------------------------|----------------|
| $\text{SO}_2$          | 0.3                          | 0.15                          |  N/A           |
| $\text{HNO}_3$         | 1.4                          | 0.5                           |  N/A           |
| $\text{HONO}$          | 1.4                          | 0.5                           |  N/A           |
| $\text{NH}_3$          | 1.4                          | 0.5                           |  N/A           |
| $\text{H}_2\text{O}_2$ | 1.4                          | 0.5                           |  N/A           |
| $\text{HCHO}$          | 0.1                          | 0.03                          |  N/A           |
| $\text{ROOH}$          | 0.05                         | 0.015                         |  N/A           |
| $\text{PM}_{2.5}$      | 1.0                          | N/A                           |  0.02          |
| $\text{PM}_{10}$       | 1.0                          | N/A                           |  0.4           |

 
## Namelist (<code>_p3d</code>) Options
 
The PALM wet deposition model can be activated from the namelist (`_p3d`) file by setting the option `chem_wet_deposition` to `.TRUE.`. The model automatically matches the available gas phase and aerosol species listed in the table above with the species defined in the active chemical mechanism in the file `chem_gasphase_mod.f90` by name.  The user can further specify the intervals in which the scavenging rates are updated via the namelist option `chem_wet_deposition_update_interval`, which is defaulted to 300 seconds.

Normally, the wet deposition model extracts information cloud layer and precipitation from the build cloud model at each update interval (defined by the `chem_wet_deposition_update_interval` option) to re-calculate the scavenging rates.  This dependence can be disabled by setting the
namelist option `chem_wet_deposition_override` to `.TRUE.`.  In this way, the vertical cloud layer and precipitation rate can be supplied directly by the user.  The location of the cloud layer is specified through the options `chem_wet_deposition_cloud_level_lower` and `chem_wet_deposition_cloud_level_uppers`, and they are defaulted to 10 and 15 respectively, denoting the vertical grid level for the lower and upper bounds of the cloud layer. Meanwhile, the precipitation rate can be specified by the option `chem_wet_deposition_rain_rate`, which carries a default of 1.0 $\text{kg}$ $\text{m}^{-2}$ $\text{s}^{-1}$ per kilogram of air.

A summary of the namelist options are listed below:

| `_p3d_` Option | Default | Function |
|----------------------------|---------|----------|
| `chem_wet_deposition` | `.FALSE.` | `.TRUE.` activates the wet deposition model |
| `chem_wet_deposition_update_interval` | 300 | Seconds between scavenging rate updates |
| `chem_wet_deposition_model_override`  | `.FALSE.` | `.TRUE.` to use user-supplied cloud and precipitation information |
| `chem_wet_deposition_cloud_level_lower` | 10 | User-supplied lower bound for the domain cloud layer |
| `chem_wet_deposition_cloud_level_upper` | 15 | User-supplied Upper bound for the domain cloud layer |
| `chem_wet_deposition_rain_rate` | 1.0 | User-supplied rate of rain precipitation for the whole domain  |

## References

- Berge, E., and Jakobsen, H.A.: A regional scale multi-layer model for the calculation of long-term transport and deposition of pollution in Europe, Tellus B, 50, 205-224, 1998.

- Laakso, L. *et al*: Ultraﬁne particle scavenging coefﬁcients calculated from 6 years ﬁeld measurements, Atmos.
Environ., 37, 3605–3613, 2003.

- Simpson, D. *et al*: The EMEP MSC-W chemical transport model – technical description, Atmospheric Chemistry and Physics, 12, 7825–7865, 2012.

- van Zanten, M.C. *et al*: Controls on precipitation and cloudiness in simulations of trade-wind cumulus as observed during RICO, Journal of Advances in Modeling Earth System 3(2), 2011.
