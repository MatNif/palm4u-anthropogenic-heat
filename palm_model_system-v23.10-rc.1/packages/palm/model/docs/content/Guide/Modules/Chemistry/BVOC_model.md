# Biogenic VOC Emission Model

---

This is the main page of the biogenic emission model (BEM) for the PALM model system. These pages describes the physical bases, structure, parameterizations, and features of BEM along with the list of pre-requisete, input data and controling parameters to run the biogenic emission model.

## Overview

The BEM for PALM computes emissions of biogenic volatile organic compounds (BVOCs) from trees and all resolved vegetation available in the model domain. The biogenic emission parameterizations in BEM are primarily based on MEGAN (Model of Emission of Gasses and Aerosols from Nature), ([Guenther et al., 1991](https://doi.org/10.1029/91jd00960); [1993](https://doi.org/10.1029/93jd00527); and [2012](https://doi.org/10.5194/gmd-5-1471-2012) ). However, certain enhancements have been made in solar radiation absorption and soil moisture impact algorithms. In order to exploit the flexibility and accuracy of the LES approach, volumetric emissions of BVOC are available for all resolved vegetation types considered. This version of BEM does not calculate emissions from unresolved (flat vegetation). The current implementation of the BEM is LOD 0 (level of detail 0).

The BVOC model implemented in PALM computes the net primary emissions of up to 30 biogenic VOCs grouped in five classes, namely isoprene, monoterpenes, sesquiterpenes, reactive oxygenated VOCs (XVOC) and other VOCs (OVOC) <a href="#T01">Table 1</a>. The net emission of each compound, F$_{i}$ (in $\mu$mole m$^{-3}$ s$^{-1}$) **from a given plant species is estimated for each vertical layer in the resolved vegetation/tree and can be described by the following equation:**

$$
      F_{i} =\sum_{j=1}^n \Bigl( \epsilon_{(i,j)} .\gamma_{i}. \Psi_{j}\Bigr)
$$

Where $\epsilon_{(i,j)}$ is the plant specific average emission potential ($\mu$mole m$^{-2}$ s$^{-1}$), of the biogenic compound class $i$ of the plant species $j$ at standard temperature of 303.15 (K) and standard photosynthetic photon flux density (PPFD) of 1000 $\mu$ mol m$^{-2}$ s$^{-1}$, $\Psi$, is the foliar density (LAD) in g dry mass m$^{-2}$) and $\gamma_{i}$ is the dimensionless emission activity factor to account for variations due to specific physiological and phenological states. Details of the calculation of $\gamma_{i}$ are given in the following section. The minimum requirement to steer BEM, is a list of emitted biogenic species. If tree types/PFTs and respective emission potentials of the emitted species are not provided by the user, BEM would read the default information from the biogenic data file. Emitted species defined by the user in the namelist, but not defined in the chosen chemistry mechanism would be ignored. The output concentrations of biogenic VOCs are in units of parts per million (ppm).

Following models are required to run the BEM,  <br/>
a) chemistry_model_mod <br/>
b) land_surface_model  <br/>
c) radiation_model_mod <br/>
d) urban_surface_model  <br/>
e) plant_canopy_model <br/>

**static_PIDS:** The BVOC model requires input vegetation data including leaf area density (LAD), tree_types, and vegetation types. The user needs to make sure that these fields exist in the static_PIDS file.


<b id="T01">Table 1</b>  : List of biogenic compounds with their respective light dependent fraction (LDF) and $\beta$ values.

<table style="width:82%;">
<colgroup>
<col style="width: 30%" />
<col style="width: 10%" />
<col style="width: 10%" />
<col style="width: 30%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Biogenic compound class</th>
<th style="text-align: left;">LDF</th>
<th style="text-align: left;">&beta;</th>
<th style="text-align: center;">Biogenic compound name</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">Isoprene</td>
<td style="text-align: center;">1.00</td>
<td style="text-align: center;">0.13</td>
<td style="text-align: left;">Isoprene</td>
</tr>
<tr class="even">
<td style="text-align: left;">Monoterpene</td>
<td style="text-align: center;">0.60</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">myrcene</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.60</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">sabinene</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">limonene</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">3-carene</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.80</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">t-&beta;-ocimen</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.60</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">&alpha;-Pinene</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.60</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">β-Pinene</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Sesquiterpenes</td>
<td style="text-align: center;">0.50</td>
<td style="text-align: center;">0.17</td>
<td style="text-align: left;">&beta;-caryophyllene</td>
</tr>
<tr class="even">
<td style="text-align: left;">Oxygenated VOCs</td>
<td style="text-align: center;">0.80</td>
<td style="text-align: center;">0.13</td>
<td style="text-align: left;">acetaldehyde</td>
</tr>
<tr class="odd">
<td style="text-align: left;"> </td>
<td style="text-align: center;">0.80</td>
<td style="text-align: center;">0.13</td>
<td style="text-align: left;">ethanol</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.80</td>
<td style="text-align: center;">0.13</td>
<td style="text-align: left;">formaldehyde</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.80</td>
<td style="text-align: center;">0.13</td>
<td style="text-align: left;">methanol</td>
</tr>
<tr class="even">
<td style="text-align: left;"> </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.13</td>
<td style="text-align: left;">acetone</td>
</tr>
<tr class="odd">
<td style="text-align: left;"> </td>
<td style="text-align: center;">0.80</td>
<td style="text-align: center;">0.13</td>
<td style="text-align: left;">farmic acid</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.80</td>
<td style="text-align: center;">0.13</td>
<td style="text-align: left;">acetic acid</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Other VOCs</td>
<td style="text-align: center;">1.00</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">232-MBO</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">methane  </td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">ethane  </td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">hydrogen cyanide</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">toluene</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">methyl bromide</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">methyl chloride</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">methyl iodide</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">dimethyle sulfide</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">ethane</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">propane  </td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">propene</td>
</tr>
<tr class="odd">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">butane</td>
</tr>
<tr class="even">
<td style="text-align: left;">  </td>
<td style="text-align: center;">0.20</td>
<td style="text-align: center;">0.10</td>
<td style="text-align: left;">benzaldehyde</td>
</tr>
</tbody>
</table>

## Emission Activity Factor ($\gamma_{i}$)

The emission activity factor $\gamma_{i}$ represents processes that cause variations in emissions due to the dependency to light and temperature, soil moisture, and general seasonality. The isoprene emission is mainly controlled by temperature and visible light within the range of 400 - 700 nm wave length ([Guenther et al., 1995](https://doi.org/10.1029/94JD02950) ), nevertheless foliage density, soil moisture, and seasonal variations also play important role ([Guenther et al., 1997](https://doi.org/10.1890/1051-0761(1997)007[0034:SASVIN]2.0.CO;2);[ 2006](https://doi.org/https://acp.copernicus.org/articles/6/3181/2006/); [and 2012](https://doi.org/10.5194/gmd-5-1471-2012) ). Monoterpene emissions from a large variety of plants originate from volatilization of stored VOCs for example from resins vessels ([Guenther et al., 1993](https://doi.org/10.1029/93jd00527)). However, monoterpene emissions from some plants are also sensitive to light <a href="#steinbrecher_1989">Steinbrecher et al., (1989)</a>, and [Guenther et al., 1995](https://doi.org/10.1029/94JD02950), for example, some Oak species (<a href="#steinbrecher_1996">Steinbrecher et al., (1996)</a>; [Seufert et al., 1997](https://doi.org/10.1016/S1352-2310(97)00334-8)). Scot pine or Norway spruce emit monoterpenes from storage pools as well as from $de$ $novo$ direct synthesis <a href="#steinbrecher_1994">(Steinbrecher et al., 1994</a>, and <a href="#steinbrecher_1999">1999)</a>. BEM for PALM considered both approaches. The calculation of emission activity in BEM is numerically described as follows:

$$
 \gamma_{i} = \left[(1 - LDF_{i}).\gamma_{(LI,i)} +  LDF{i} .\gamma_{(P)}.\gamma_{(T)}.\gamma_{(SM)}. \gamma_{(SN)} \right]
$$

Where $\gamma_{(P)}$, and $\gamma_{(T)}$ are dimensionless emission activity factors for light and temperature correction, respectively for the newly synthesised monoterpenes and isoprene, and $\gamma_{(LI)}$ is the light independent emission factor that follow the monoterpene exponential temperature response described by [Guenther et al., (1993)](https://doi.org/10.1029/93jd00527), $\gamma_{(SN)}$, is the emission activity factor for the seasonal influence, and $\gamma_{(SM)}$ is the soil moisture dependence (calculated in m$^{3}$ m$^{-3}$). Emission of each biogenic compound includes a light dependent fraction (LDF) <a href="#T01">Table 1</a> with the remaining light independent fraction (1-LDF) that is not influenced by light. <a href="#T01">Table 1</a>:

**a) Light and temperature dependence of biogenic emissions**
The light dependent emission activity is calculated following the isoprene response to temperature as described by [Guenther et al., (1993)](https://doi.org/10.1029/93jd00527):

$$
 \gamma_{(P)}  = \left[\frac{\alpha\,C_{L1}L } {\sqrt{1 + \alpha^2 L^2} } \right]
$$

where L is photosynthetically active photon flux density (PPFD) in $\mu$mol photons m$^{-2}$ s$^{-1}$, and $C_{L1}$ (=
1.066), and $\alpha$ (= 0.0027) are empirical coefficients [Guenther et al., (1993)](https://doi.org/10.1029/93jd00527). PPFD (in $\mu$mol photons m$^{-2}$ s$^{-1}$) indicates photosynthetically active radiation (PAR) within the visible range (400-700nm) absorbed by plants. PPFD can be estimated by multiplying total shortwave radiation (direct + diffused) with an empirical factor of
2.02   ([Reis et al., 2020](https://doi.org/10.31062/agrom.v27i2.26527)). Three methods are available to calculate the total shortwave radiation incident on a grid-cell. These are described in the namelist parameter list under ebio_rad_method (<a href="#T02">Table 2</a>).

The temperature dependence of biogenic emissions for BEM is given by:

$$
    \gamma_{(T)} = \left[\frac{exp \left(\frac{C_{T1}(T - T_s)}{RT_s T}\right)}{C_{T3}+exp \left(\frac{C_{T2}(T - T_M)}{RT_s T} \right)} \right]
$$

where T is the leaf temperature in Kelvin. Following [Lindfors and Laurila (2000)](http://www.borenv.net/BER/archive/pdfs/ber5/ber5-095.pdf), and others, in this version of the BVOC emission  model, we set the leaf temperature equal to the ambient temperature of the respective plant canopy layer. Parameter 'R' is the gas constant (=8.314 J K$^{-1}$ mol$^{-1}$), $C_{T1}$ (=95000 J mol$^{-1}$), $C_{T2}$ (=230000 J mol$^{-1}$), $C_{T3}$ (= 0.961), and $T_M$ (= 314 K), are empirical coefficients.

Short term variations in monoterpene emission rates from storage change with the vapour pressure of essential oils that in turn is determined by leaf temperature ([Geron et al.,  1994]( https://doi.org/https://doi.org/10.1029/94JD00246)). For volatilization-controlled emission,[Guenther et al., (1993)](https://doi.org/10.1029/93jd00527) recommends the following emission response:

$$
 \gamma_{(LI,i)} = exp(\beta_i[T -  T_s])
$$

where $\beta$ is an empirical coefficient specific to the biogenic compound class (see <a href="#T01">Table 1</a>) and T${_s}$ is the leaf temperature in standard conditions (=303 K)

**b) Soil moisture ($\gamma_{SM}$)**
Biogenic emissions, in particular isoprene, are directly and indirectly influenced by soil moisture. This generally
results in emission reduction as soon as a soil moisture threshold is crossed  [(Guenther et al., 2006)](https://doi.org/https://acp.copernicus.org/articles/6/3181/2006/). Accordingly, an emission activity factor that describes this soil moisture dependence $\gamma_{SM}$ has been included. Two different approaches to calculate effective soil moisture have been implemented: a) a simplified $'bulk'$ parameterization which is also used in the Guenther model [(Guenther et al., 2012)](https://doi.org/10.5194/gmd-5-1471-2012), and b) a more detailed approach that considers the distribution of fine roots and soil properties with depth (the $'weighted'$ parameterization). Currently, the soil moisture dependence algorithm is explicitly calculated for isoprene only; for all other compounds $\gamma_{SM}$ = 1. The $'bulk'$ dependence is calculated
as follows:

\begin{aligned}
& \gamma_{(SM)} = 1    \hspace{7.2cm} \theta > \theta_{1}\\
& \gamma_{(SM)} = (\theta - \theta_{w})/ \Delta \theta_{1} \hspace{5.3cm} \theta_{w} < \theta < \theta_{1}\\
&\gamma_{(SM)} =  0  \hspace{7.10cm}     \theta <  \theta_{w}
\end{aligned}

Where $\theta$ is overall soil moisture (volumetric water content, m$^3$ m$^{-3}$) averaged over the whole soil profile, $\theta_{w}$ is the wilting point (soil moisture level below which plants cannot extract water from soil, m$^3$m$^{-3}$), $\Delta \theta_{1}$ is an empirically set parameter (0.04), so that the water holding capacity ($\theta_{1}$) is assumed as $\theta_{w} + \Delta \theta_{1}$. The wilting point depends on soil type. It is lower in sandy than in a loamy soil and ranges typically between 10 and 20 percent of soil volume (see e.g. [Qi et al., 2018](https://doi.org/10.1002/hyp.11452)). [Guenther et al., (2006)](https://doi.org/https://acp.copernicus.org/articles/6/3181/2006/), used the global wilting point dataset from [Chen and Dudhia (2001)](https://journals.ametsoc.org/view/journals/mwre/129/4/1520-0493_2001_129_0569_caalsh_2.0.co_2.xml?tab_body=fulltext-display)  while the soil moisture algorithm in [Guenther et al., 2012](https://doi.org/10.5194/gmd-5-1471-2012) is based on [Pegoraro et al., (2004)](https://doi.org/10.1016/j.atmosenv.2004.07.028). The BVOC model accounts for soil moisture fro ECMWF-IFS described on PALM web pages under [soil type](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#soil_type)

In the 'weighted' approach, soil moisture influence accounts for the field capacity and wilting point as well as root distribution along the rooted soil profile to define the drought stress factor for individual trees and/or plant canopy. In this procedure $\gamma_{(DS)}$ is calculated as follows:

$$
\gamma_{(DS)} =\sum_{k=0}^n \Bigl( soilf_{(k)} + ( rootf_{(k)} * 0.5 * rwc_{(k)} ) \Bigr)
$$

Where $\gamma_{(DS)}$ is the drought stress factor, soilf is the soil fraction calculated from [deep soil layers](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#dz_soil), $rootf$ is the [root fraction](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#root_fraction), and $rwc$ is the relative water content which is defined as:

$$
rwc_{(k)} = ( soilm_{(k)} - qw_{(k)} ) / ( qc_{(k)} - qw_{(k)} )
$$

Where $soilm_{(k)}$ is the soil moisture content at soil layer $k$, and $qw_{(k)}$ and $qc_{(k)}$ are [wilting point](http://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#soil_type), and [field capacity](http://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#soil_type), respectively at soil layer $k$. Both procedures can be selected in the BEM.

**c) Seasonal Influences ($\gamma_{SN}$)**
Plants emit different amounts of isoprene even under the same light and temperature conditions depending on the day of the year (or the activity stage of the foliage, respectively)  ([Ohta, 1986](https://doi.org/https://www.jstage.jst.go.jp/article/geochemj1966/19/5/19_5_269/_article/-char/ja/); [Monson et al., 1994](https://doi.org/10.1007/bf00627738) ). Isoprene emissions follow a seasonal pattern of negligible emissions in the winter period followed by a rapid rise, increasing to a growing season maximum and a decrease towards winter ([Monson et al., 1994](https://doi.org/10.1007/bf00627738); [Goldstein, 1994](https://www.proquest.com/openview/54cf57b04e344719e4fcf0c404e54152/1?pq-origsite=gscholar&cbl=18750&diss=y);[Guenther et al., 1997](https://doi.org/10.1890/1051-0761(1997)007[0034:SASVIN]2.0.CO;2);  [Kempf et al., , 1996](https://doi.org/10.1016/1352-2310(95)00462-9). [Staudt et al., 2000]( https://doi.org/10.1023/A:1006233010748) ), proposed a symmetric equation describing a Gaussian (bell shaped) response in order to modify light-dependent monoterpene as well as isoprene emissions. [Keenan et al., (2009)](https://doi.org/10.1029/2009jd011904) compared different shapes of seasonality functions and concluded that an asymmetric function better adheres to the data and is recommended for simulation of seasonal variation of isoprenoid emission [Grote et al., (2013)](https://link.springer.com/chapter/10.1007/978-94-007-6606-8_12) . [Keenan et al., (2009)](https://doi.org/10.1029/2009jd011904) directly adjusted the emission rate $F_i$, rather than creating a multiplier of $\epsilon_i$. We therefore used a modified version of the asymmetric function proposed by [Keenan et al., (2009)](https://doi.org/10.1029/2009jd011904). The general behaviour of the seasonal dependence of biogenic emissions is thus numericallyd escribed as follows:

$$
\gamma_{SN} = exp[ - (D - D_{max})/\tau)]^2
$$

Where $D$ is the day of the year, $D_{max}$(=200) represents the day on which the emission capacity reaches its maximum and $\tau$ (=100) is the breadth (kurtosis) of the seasonal amplitude (in days). Both parameters are adjustable through the parameter namelists. It is to be noted that the mathematical form of $\gamma_{SN}$ does not represent changes of leaf abundance during the season that also affect $\epsilon_i$ or $Fi$, because this influence is already reflected in the determination of foliar density. Assuming further that storage pools are not depleted during the year, the seasonal variation factor only modifies the light dependent emissions. For emissions from storage, thus $\gamma_{SN}$ = 1.


## INPUT PARAMETERS

<a href="#T02">Table 2</a> provides list of all namelist parameters that may be required to steer the BEM. All namelist input parameters should be defined in the chemistry namelist &chemistry_parameters.

All BVOC output quantities follow the same convention as applicable to other chemical species and aerosols i.e. $kc\_$ followed by the species name in the respective mechanism.


###  `chemistry_parameters` namelist

<b id="T02">Table 2</b>: List of input parameters for the biogenic emission model.

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
<td style="text-align: left;">ebio_dt</td>
<td style="text-align: center;">R</td>
<td style="text-align: center;">'0.0'</td>
<td style="text-align: left;">Time-step of biogenic emission model, that is how often biogenic emission would be updated. The default value is 0. If user does not define ebio_dt, the BEM will use models physical time-step. In case user defines ebio_dt but it is less than model's physical time-step then model will issue a warning and make ebio_dt = dt_timestep.</td>
</tr>
<tr class="even">
<td style="text-align: left;">ebio_ef_pft</td>
<td style="text-align: center;">R</td>
<td style="text-align: center;">9999999.9</td>
<td style="text-align: left;"> emission potentials for each 
<a href="https://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#vegetation_type" target=_blank>vegetation type (PFT)</a> </td>
</tr>
<tr class="odd">
<td style="text-align: left;">ebio_ef_tree</td>
<td style="text-align: center;">R</td>
<td style="text-align: center;">9999999.9</td>
<td style="text-align: left;"> emission potentials for each <a href="https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids/palm_csd#Singletrees" target=_blank>tree type (single trees)</a> </td>
</tr>
<tr class="even">
<td style="text-align: left;">ebio_emis_name</td>
<td style="text-align: center;">C</td>
<td style="text-align: center;">'novalue'</td>
<td style="text-align: left;">  Names of biogenic VOCs. The species name should be available in both BVOC <a href="#T02">Table 2</a> as well as in the chemical mechanism. User should follow the naming convention in th given chemistry mechanism. At the moment only CBM4 and SMOG chemistry  mechanisms have been implemented in the biogenic model.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ebio_max_emis_day</td>
<td style="text-align: center;">I</td>
<td style="text-align: center;">200</td>
<td style="text-align: left;"> The day when the emissions are expected to be the maximum, The default is 200 which is 19<sup>th</sup> or 20<sup>th</sup> of July.</td>
</tr>
<tr class="even">
<td style="text-align: left;">ebio_pft</td>
<td style="text-align: center;">I</td>
<td style="text-align: center;">-127</td>
<td style="text-align: left;">Vegetation types (plant functional types). PALM utilizes ECMWF-IFS classification for <a href="https://palm.muk.uni-hannover.de/trac/wiki/doc/app/land_surface_parameters#vegetation_type">vegetation types</a> </td>
</tr>
<tr class="odd">
<td style="text-align: left;">ebio_ppfd_factor</td>
<td style="text-align: center;">R</td>
<td style="text-align: center;">202</td>
<td style="text-align: left;"> Photosynthetic photo flux density factor. The default value is the conversion factor of solar radiation (direct and diffused) incident on a grid box from Wm<sup>-2</sup> to mol m<sup>-2</sup> s<sup>-1</sup></td>
</tr>
<tr class="even">
<td style="text-align: left;">ebio_rad_method</td>
<td style="text-align: center;">I</td>
<td style="text-align: center;">0</td>
<td style="text-align: left;">Radiation method for biogenic model defines as how radiation would be used to calculate bvoc emissions. Two methods have been planned to estimate the flux density received by each PCB. The first method(= 0), (default) is based on calculating the average shortwave radiation (both direct and diffuse components) at each PCB, considering all radiative transfer processes, including shade and reflections. The RTM calculates this flux density using the absorbed radiation at each PCB and, given that albedo is set to zero (see <a href="https://doi.org/10.5194/gmd-14-3095-2021">(Pavel et al., 2021)</a>. The second planned method (= 1), is based on calculating the shade ratio of each PCB and the radiation received for the unshaded part. This method takes into account that the radiation reduces from the top face toward the inner centre of the PCB because the leaves on the top surface receive full radiation and leaves inside the grid cell receive less radiation due to being in the shade of the leaves above. This method is not yet implemented. </td>
</tr>
<tr class="odd">
<td style="text-align: left;">ebio_soilm_method</td>
<td style="text-align: center;">C</td>
<td style="text-align: center;">'bulk'</td>
<td style="text-align: left;">The BEM offers the choice of algorithm in calculation of soil moisture correction factor. The default option bulk is the simple method same as implemented in MEGAN <a href="https://doi.org/10.5194/gmd-5-1471-2012">(Guenther et al., 2012) </a>, whereas weighted' is more complex and needs detailed data of soil, wilting point and field capacity</td>
</tr>
<tr class="even">
<td style="text-align: left;">ebio_tree</td>
<td style="text-align: center;">I</td>
<td style="text-align: center;">-127</td>
<td style="text-align: left;">Tree types in the simulation domain. The current PALM database has 86
<a href="https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids/palm_csd#Singletrees">tree types </a>. The user need to choose tree types from this table.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">emis_biogenic</td>
<td style="text-align: center;">L</td>
<td style="text-align: center;">.FALSE.</td>
<td style="text-align: left;">Parameter to switch biogenic emission on/off.</td>
</tr>
<tr class="even">
<td style="text-align: left;">emis_biogenic_lod</td>
<td style="text-align: center;">I</td>
<td style="text-align: center;">0</td>
<td style="text-align: left;">Level of detail i.e. biogenic emission mode. The BEM is the default LOD = 0 model and currently the
only implemented mode. LOD-1 is MEGAN 2.3 and LOD-2 is the advanced BEM (ABEM). Both LOD-1 and LOD-2 are not implemented as yet.</td>
</tr>
</tbody>
</table>

## EMISSION INPUT

The BEM for PALM does not require any emission input file for LOD0. All emission input is provided either through the chemistry namelist or in the chem_emis_biogenic_data_mod file. The BEM LOD1 is MEGAN 2.3 coupled to the PALM model. The BEM-LOD 1 is not yet implemented, however, additional input file(s) will be required for the BEM-LOD1. The BEM LOD2 is the advanced biogenic emission model (ABEM), also not implemented as yet. ABEM will also require additonal input file(s).

### EXAMPLE SETUP

Files to carry out example runs with biogenic emissions turned on.

### Example namelist file

```

&chemistry_parameters

chem_gasphase_on          = .TRUE.,
emis_generic              = .TRUE.,
chem_mechanism            = 'cbm4', 
deposition_dry            = .TRUE.,
photolysis_scheme         = 'simple',

emis_biogenic             = .TRUE.,       ! activate bvoc model,  default is .FALSE.
emis_biogenic_lod         = 0,            ! 0 = BEM (default), 1= MEGAN2.3, 2= ABEM
                                          ! only default is available.
ebio_tree                 = 2, 59,60,61,  ! different single treetypes
ebio_pft                  = 0, -7,-9      ! 0 is default pft/tree type, -7 and -9 are
                                          ! the vegetation type 7 and 9.
ebio_soilm_method         = 'weighted',   !1. 'bulk'(default), 2. 'weighted'
ebio_ppfd_factor          = 200,          ! day of the maximum emission of the bvoc.
ebio_dt                   = 0.5,          ! time step (in seconds) for the bvoc model  
ebio_rad_method           = 0,            ! 0:default, 1: absorption, 2:Fractional, 
                                          ! only default is available
    
ebio_emis_name            = 'ISOP','HCHO','CO','ETH','TOL', ! As defined in the respective chem mech 
                                                            ! in this case it is cbm4, 5 categories: 
                                                            ! isop, mterp,sqtrp, xvoc,ovoc species
                                                            ! that are not definedin the chem mechanism,
                                                            ! willl be ignored.
ebio_ef_pft(1,:)          = 0.00011, 0.00062, 0.00023,      ! Emission factor of vegetation types/PFTs
ebio_ef_pft(2,:)          = 0.00021, 0.00072, 0.00033,      ! of the form ( biogenic_spcs, pft )
                                                            ! in umol m-2 sec-1
ebio_ef_pft(3,:)          = 0.00031, 0.00082, 0.00043,
ebio_ef_pft(4,:)          = 0.00041, 0.00092, 0.00053,
ebio_ef_pft(5,:)          = 0.00051, 0.00012, 0.00063,

ebio_ef_tree(1,:)         = 0.00035,  0.00030,  0.00020,  0.00010,  ! Emission potentials/factors of trees 
ebio_ef_tree(2,:)         = 0.00041,  0.00035,  0.00025,  0.00075,  ! trees of the form (biogenic_spcs,  
ebio_ef_tree(3,:)         = 0.00045,  0.00039,  0.00055,  0.00080,  ! tree types) in umol m-2 sec-1
ebio_ef_tree(4,:)         = 0.00050,  0.00036,  0.00060,  0.00085,
ebio_ef_tree(5,:)         = 0.00015,  0.00040,  0.00065,  0.00095,

/

```
## CHEMISTRY MECHANISMS

Currently only two of the default chemistry mechanisms (CBM4 and SMOG) have some BVOCs and therefore can be used with the BEM. More chemistry mechanisms such as RADM and CBMZ are planned to implement in the near future. List of all chemistry mechanisms currently available in PALM-4U are linked to documentation. Detailed instructions as how to create a new gas-phase mechanism is available on [chemistry mechanisms](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/chemmech) page.

## LIMITATIONS

The current version of the BEM can only compute emissions from resolved (high) vegetation with leaf area density (LAD) available. All flat vegetation without LAD is ignored. The BEM is able to compute emissions of only those compounds that are available in the chosen chemical mechanism. Any VOC that is in the compound list but not in the mechanism or vice-a-versa, would be ignored.

The user is responsible to provide vegetation information in the format as described on [static driver](https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids/static) page. Compound class specific model parameters for the calculation of light and temperature activity factor; adopted from [(Henrot et al., 2017)](https://doi.org/10.5194/gmd-10-903-2017)


## BIBLIOGRAPHY

Chen, F., Dudhia, J., 2001. Coupling and advanced land surface-hydrology model with the Penn State-NCAR MM5 modeling system. Part I: Model implementation and sensitivity. Mon. Weather Rev. 129, 569–585. https://doi.org/10.1175/1520-0493(2001)129<0569:CAALSH>2.0.CO;2

Forkel, R., Klemm, O., Graus, M., Rappenglück, B., Stockwell, W.R., Grabmer, W., Held, A., Hansel, A., Steinbrecher, R., 2006. Trace gas exchange and gas phase chemistry in a Norway spruce forest: A study with a coupled 1-dimensional canopy atmospheric chemistry emission model. Atmos. Environ. 40, 28–42. https://doi.org/10.1016/j.atmosenv.2005.11.070


Geron, C.D., Guenther, A.B., E., P.T., 1994. An improved model for estimating emissions of volatile organic compounds from forests in the eastern United States. Jounal Geophys. Res. 99, 12773–12791. https://doi.org/https://doi.org/10.1029/94JD00246


Goldstein, A., 1994. Non-methane hydrocarbons above a midlatitude forest: Biogenic emissions and seasonal concentration variations.Harvard University. https://www.proquest.com/openview/54cf57b04e344719e4fcf0c404e54152/1?pq-origsite=gscholar&cbl=18750&diss=y

Grote, R., Monson, R.K. and Niinemets, Ü., 2013. Leaf-level models of constitutive and stress-driven volatile organic compound emissions. In Biology, controls and models of tree volatile organic compound emissions (pp. 315-355). Springer, Dordrecht.https://link.springer.com/chapter/10.1007/978-94-007-6606-8_12

Guenther, A., 1995. A global model of natural volatile organic compound emissions. J. Geophys. Res. 100, 8873–8892. https://doi.org/10.1029/94JD02950

Guenther, A., 1997. Seasonal and spatial variations in natural volatile organic compound emissions. Ecol. Appl. 7, 34–45. https://doi.org/10.1890/1051-0761(1997)007[0034:SASVIN]2.0.CO;2

Guenther, A.B., Karl, T., Harley, P., Wiedinmyer, C., Palmer, P.I., Geron, C., 2006. Estimates of global terrestrial isoprene emissions using MEGAN (Model ofEmissions of Gases and Aerosols from Nature). Atmos. Chem. Phys, 6, 3181–3210. https://doi.org/https://acp.copernicus.org/articles/6/3181/2006/


Guenther, A.B., Karl, T., Harley, P., Wiedinmyer, C., Palmer, P.I., Geron, C., 2006. Estimates of global terrestrial isoprene emissions using MEGAN (Model ofEmissions of Gases and Aerosols from Nature). Atmos. Chem. Phys, 6, 3181–3210. https://doi.org/https://acp.copernicus.org/articles/6/3181/2006/

Guenther, A.B., Jiang, X., Heald, C.L., Sakulyanontvittaya, T., Duhl, T., Emmons, L.K., Wang, X., 2012. The model of emissions of gases and aerosols from nature version 2.1 (MEGAN2.1): An extended and updated framework for modeling biogenic emissions. Geosci. Model Dev. 5, 1471–1492. https://doi.org/10.5194/gmd-5-1471-2012

Guenther, A.B., Monson, R.K., Fall, R., 1991. Isoprene and monoterpene emission rate variability: Observations with eucalyptus and emission rate algorithm development. J. Geophys. Res. 96, 10799. https://doi.org/10.1029/91jd00960

Guenther, A.B., Zimmerman, P.R., Harley, P.C., Monson, R.K., Fall, R., 1993. Isoprene and monoterpene emission rate variability: model evaluations and sensitivity analyses. J. Geophys. Res. 98. https://doi.org/10.1029/93jd00527

Henrot, A.J., Stanelle, T., Schröder, S., Siegenthaler, C., Taraborrelli, D., Schultz, M.G., 2017. Implementation of the MEGAN (v2.1) biogenic emission model in the ECHAM6-HAMMOZ chemistry climate model. Geosci. Model Dev. 10, 903–926. https://doi.org/10.5194/gmd-10-903-2017

Keenan, T., Niinemets, Ü., Sabate, S., Gracia, C., Peñuelas, J., 2009. Seasonality of monoterpene emission potentials in Quercus ilex and Pinus pinea : Implications for regional VOC emissions modeling . J. Geophys. Res. 114, 1–11. https://doi.org/10.1029/2009jd011904

Kempf, K., Allwine, E., Westberg, H., Claiborn, C., Lamb, B., 1996. Hydrocarbon emissions from spruce species using environmental chamber and branch enclosure methods. Atmos. Environ. 30, 1381–1389. https://doi.org/10.1016/1352-2310(95)00462-9

Lindfors, V., Laurila, T., 2000. Biogenic volatile organic compound (VOC) emissions from forests in Finland. Boreal Environ. Res. 5, 95–113. http://www.borenv.net/BER/archive/pdfs/ber5/ber5-095.pdf

Monson, R., Harley, P., Litvak, M., Oecologia, M.W.-, 1994, undefined, 1994. Environmental and developmental controls over the seasonal pattern of isoprene emission from aspen leaves. Springer 99, 260–270. https://doi.org/10.1007/bf00627738

Ohta, K., 1986. Diurnal and seasonal variations in isoprene emission from live oak. Geochem. J. 19, 269–274. https://doi.org/https://www.jstage.jst.go.jp/article/geochemj1966/19/5/19_5_269/_article/-char/ja/

Krč, P., Resler, J., Sühring, M., Schubert, S., Salim, M. H., & Fuka, V., 2021. Radiative Transfer Model 3.0 integrated into the PALM model system 6.0. Geoscientific Model Development, 14(5), 3095-3120. https://doi.org/10.5194/gmd-14-3095-2021

Pegoraro, E., Rey, A., Greenberg, J., … P.H.-A., 2004, 2004. Effect of drought on isoprene emission rates from leaves of Quercus virginiana Mill. Elsevier. https://doi.org/10.1016/j.atmosenv.2004.07.028

Qi, J., Markewitz, D., Radcliffe, D., 2018. Modelling the effect of changing precipitation inputs on deep soil water utilization. Hydrol. Process. 32, 672–686. https://doi.org/10.1002/hyp.11452

Reis, M.G. dos, Ribeiro, A., 2020. Conversion factors and general equations applied in agricultural and forest meteorology. Agrometeoros 27, 227–258. https://doi.org/10.31062/agrom.v27i2.26527

Seufert, G., Bartzis, J., Bomboi, T., Ciccioli, P., Cieslik, S., Dlugi, R., Foster, P., Hewitt, C.N., Kesselmeier, J., Kotzias, D. and Lenz, R., 1997. An overview of the Castelporziano experiments. Atmospheric Environment, 31, 5-17. https://doi.org/10.1016/S1352-2310(97)00334-8

Staudt, M., Bertin, N., Frenzel, B., Seufert, G., 2000. Seasonal variation in amount and composition of monoterpenes emitted by young Pinus pinea trees - Implications for emission modeling. J. Atmos. Chem. 35, 77–99. https://doi.org/10.1023/A:1006233010748

<p id="steinbrecher_1989">Steinbrecher, R. (1989). Gehalt und Emission von Monoterpenen in oberirdischen Organen von Picea abies. Ph.D. thesis, Tech. Univ. Munchen, Germany.

<p id="steinbrecher_1994">Steinbrecher, R. (1994). Emission of vocs from selected european ecosystems: the state of the art. In The proceedings of EUROTRAC symposium, volume 94, pages 448–454.

<p id="steinbrecher_1996">Steinbrecher, R. and Hauff, K. (1996). Isoprene and monoterpene emission from mediterranean oaks. In The proceedings of EUROTRAC symposium, volume 96, pages 229–233.

<p id="steinbrecher_1999">Steinbrecher, R., Hauff, K., Hakola, H., and Rossler, J. (1999). A revised parameterisation for emission modelling of isoprenoids for boreal plants. Air Pollution research report, 70(3):2944.</p>
