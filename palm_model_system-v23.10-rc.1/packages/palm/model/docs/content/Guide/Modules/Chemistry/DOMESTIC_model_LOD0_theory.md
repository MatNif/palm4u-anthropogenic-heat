# Theory for Domestic Heating Emissions Parametrization

The theoretical foundations for parameteric modelling of domestic heating emissions are derived from the works of Baumbach et al. (2010) and Struschka and Li (2019), which are based on a direct relationship between emissions and energy consumption. Emissions are calculated using the so-called *emission factors*, which are species and technology dependent, while energy consumption are functions of the size, geometry, age, and function of the individual buildings. These forms the two aspects of the discussion below.  Typically, buildings with a footprint of less than 10 m<sup>2</sup> and an effective height of less than 3 m are not considered in the calculation. 

As values of all required parameters are either given as default or provided as documented in the [domestic model overview](./DOMESTIC_model.md), this document will focus on the mathematical foundations for the parametrized model. The reader is also encouraged to refer to the works listed in the references section below for further detail.

## Daily Energy Consumption

The daily energy consumption of a building ($E_B$) at any given time can be expressed as a ratio between its annual aggregate and the number of degrees to be heated to maintain its target indoor temperature:

$\displaystyle E_B\left(t\right) = \left(\frac{E_B}{\Delta{T}}\right)_A \Delta{T}\left(t\right){\quad\quad}$ **(1)**,

where<br>
${\quad}\left.E_B\right|_A$ is the annual energy consumption of building $B$,<br>
${\quad}\left.\Delta{T}\right|_A$ is annually accumulated temperature deficit, also known as *heating degrees*, and<br>
${\quad}\Delta{T}\left(t\right)$ is the current temperature deficit, known as *heating degree days*.

The value of $\Delta{T}_A$ is provided as user input. The calculation of the remaining terms will be further presented below.

### Annual Energy Consumption

The annual energy consumption takes into account the volume, compactness, and energy demand of the building:

$\displaystyle \left.E_B\right|_A = \left.\kappa_\beta\right|_A\Phi_\beta{V_B}{\quad\quad}$ **(2)**, 

where<br>
${\quad}\Phi_\beta$ is the compactness factor of the building type $\beta$ belonging to building $B$,<br>
${\quad}\left.\kappa_\beta\right|_A$ is the annual energy demand of the building type $\beta$ belonging to building $B$ per unit footprint area, and<br>
${\quad}V_B$ is the volume of the building $B$.<br>

The compactness factor, in unit of 1/m, is a density indicator of the building. On the other hand, the annual energy demand is the amount of energy to be consumed annually per unit area of the building. Tabulated values of these two quantities are available for each building type ($\beta$) and are used as default inputs.

### Temperature Deficit (Heating Degree Days)

The temperature deficit, or heating degree days, is calculated by subtracting the user-defined base temperature from the outdoor ambient temperature:

$\displaystyle \Delta{T}\left(t\right) = \max\left\{\ 0, \left[T_0 - T_\infty(t)\right]\ \right\}{\quad\quad}$ **(3)**,

where<br>
$T_0$ is the base temperature, and<br>
$T_\infty\left(t\right)$ is the ambient temperature.

When the ambient temperature ($T_\infty$) is greater than the base temperature ($T_0$), there will be no temperature deficit (i.e., $\Delta{T}=0$) and it is assumed that no heating will be required.

The ambient temperature is, in turn, the mean temperature over the volume for region of interest ($V$) up to a user-defined sampling height above ground level ($\eta$):

$\displaystyle T_\infty\left(t\right) = \frac{1}{V}\iiint_\eta{T}dV{\quad\quad}$ **(4)**.

Ideally, the sampling height ($\eta$) should be chosen so that the temperature is representative of that in the urban canopy. Because of this, to reduce computing resources, $\Delta{T}_D$ is only calculated at fixed time intervals to minimize computing resources.

## Emissions Source

The emission of the individual species ($\epsilon_B^k$) can then be calculated based on the building energy consumption $E_B$:

$\displaystyle \epsilon_B^k = E_B\psi^k{\quad\quad}$ **(5)**,

where $\psi$ is the emission factor for pollutant species $k$.

The emission factors $\psi$ are tabulated on either a molar or mass basis for the energy consumed, and are expressed in mol/TJ for gas phase species and kg/TJ for particle species.

## Nomenclature

*Note: Input quantities to the model are indicated in parentheses.*

| Symbol         | Unit               | Description                                     |
|----------------|--------------------|-------------------------------------------------|
| $A$            | - -                | Subscript denoting annual aggregate             |
| $B$            | - -                | Subscript denoting building                     |
| $E$            | J                  | Energy consumption                              |
| $k$            | - -                | Superscript denoting emission species           |
| $T_0$          | K                  | Base temperature (input)                        |
| $T_\infty$     | K                  | Ambient temperature                             |
| $t$            | s                  | Time                                            |
| $V$            | m<sup>3</sup>      | Volume                                          |
| $\beta$        | - -                | Subscript denoting building type                |
| $\Delta{T}$    | K                  | Temperature deficit (heating degree days)       |
| $\Delta{T}_A$  | K                  | Annual heating degrees (input)                  |
| $\epsilon$     | mol/s or kg/s      | Emission                                        |
| $\kappa$       | J/m<sup>2</sup>    | Energy demand per unit footprint area (input)   |
| $\Phi$         | 1/m                | Compactness factor (input)                      |
| $\eta$         | m                  | Sampling height for ambient temperature (input) |

## References

- Baumbach, G., Struschka, M., Juschka, W., Carrasco, M., Ang, K.B., Hu, L. (2010) Modellrechnungen zu den Immissions-belastungen bei einer verstärkten  Verfeuerung von Biomasse in Feuerungsanlagen der 1. BImSchV. Umweltbundesamt (Bessau-Roßlau), ISSN 1862-4804.

- Struschka, M., Li, L. (2019) Temperaturabhängige zeitliche Disaggregation von Emissionen aus Feuerungsanlagen der Haushalte und Industrie für Berlin im Rahmen des MOSAIK-Projektes. Universtiät Stuttgart.
