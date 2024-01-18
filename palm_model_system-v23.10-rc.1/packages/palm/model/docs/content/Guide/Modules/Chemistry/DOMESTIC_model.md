# Domestic Heating Emission Mode

The domestic heating emission mode can be activated by setting the value of the `_p3d` option `emis_domestic_heating` to `.TRUE.`.  The domestic heating mode is currently available in LOD 0 and LOD 2.  The user can select the specific LOD for the emission mode by using the `_p3d` option `emis_domestic_lod`.  The default is `emis_domestic_lod = 2`.

## Specifications under LOD 0

Emissions are parameterized under LOD 0, where the user can specify its characteristics in the `_p3d` file as well as in the `_static` file. The physical and type of all buildings, as well as the locations of the corresponding stack (chimney) will be extracted from the `_static` file during initialization. Meanwhile, in the `_p3d` file, the user can define the heating characteristics of each building type, emission species, corresponding emission factors, as well as an overall base temperature to represent energy exchange between the buildings and its outdoor surroundings.

The emissions source will be calculated once at the start of each model run, and thereafter at a user-defined update interval (in the `_p3d` file). During each update, the deficit in base and mean ambient temperatures will first be determined. Should a deficit exists (i.e., the ambient temperature is lower than the base temperature), a non-zero rate of energy consumption of each building will be further calculated based on the parametric inputs.  The emissions for each valid species (i.e., species defined in the `_p3d` file that are present in the active chemical mechanism) are calculated using corresponding emission factors.  The emissions are released in to the computational domain above each stack location as volumetric source terms. 
 
The theoretical foundations behind the parametrization can be found in the link below:

[Theory for domestic heating emissions parameterization](DOMESTIC_model_LOD0_theory.md)

### Required Variables the `_static` File

The following variables in the `_static` file must be present and defined for LOD 0:

* `buildling_2d` - Containing the footprint of the buildings and heights at each cell location.
* `building_type` - The class and function of the buildings
* `stack_building_volume` - The volume of the building reported at the location of the stacks (chimneys)

All variables are in two dimensions presenting the horizontal extent of the model interest in cell (j,i) coordinates.

Further, six building types are supported to represent their heating characteristics:

| Type | Function    | Class       | 
|------|-------------|-------------|
| 1    | Residential | Before 1950 |
| 2    | Residential | 1951 - 2000 |
| 3    | Residential | After 2000  | 
| 4    | Commercial  | Before 1950 |
| 5    | Commercial  | 1951 - 2000 |
| 6    | Commercial  | After 2000  |

The building type is defined by its *compactness factor* and *annual energy demand*; they will be explained in detail below. 
 
### Options for the `_p3d` File

All `_p3d` options on domestic heating emissions carry the prefix `emis_domestic_`.  The following options can be defined in the `_p3d` file. Options without default values are mandatory:

| `emis_domestic_` ...        | Unit                    |Default    | Description |
|-----------------------------|-------------------------|-----------|-------------|
| `update_interval`           | s                       | 300       | Interval between emission source updates |
| `sampling_k`                | None                    | 10        | Number of vertical layers above ground for ambient temperature calculation |
| `heating_degree`            | K p.a.                  | 2100      | Annual heating degree |
| `base_temperature`          | $^\circ$C               | 15        | Base temperature for determining temperature deficit |
| `compact_factors`           |  1/m                    | See below | Compactness factors for each building type |
| `energy_demands`            | kWh/m<sup>2</sup> p.a. | See below | Annual energy demand for each building type |
| `species_names`             | N/A                     | None      | Name of emission species |
| `species_emission_factors ` | mol/TJ or kg/TJ         | None      | Emission factor for each species defined |

`emis_domestic_update_interval`: Emission sources at each stack location are updated at the start of the model run, and afterwards at every user-defined update interval. The default is 300 seconds, but any positive value can be specified.

`emis_domestic_sampling_k`: The ambient temperature for each computational domain are calculated during each update by taking the average temperature of all the cells in the domain up to a certain vertical layers above ground, defined by the option `emis_domestic_sample_k`.  The default value is 10.

`emis_domestic_heating_degree`: The amount of heating required on an annual basis are inferred by the so-called *annual heating degree*, that is, the cumulative temperature, in degrees, to be heated above ambient temperature to the base temperature throughout the year. A value of 2100 K is presented as default.

`emis_domestic_base_temperature`: The mean domestic indoor temperature to be maintained by heating. The default is 15 $^\circ$C.

`emis_domestic_compact_factors` and `emis_domestic_energy_demands`: Compactness factor is an indicator of building density, while the annual specific energy demand estimates the energy requirement for the building per unit footprint area. These building characteristics are defined for each building type (1-6, see section on `_static` file). The table below shows the default values used, which can be individually overwritten by the user in the `_p3d` file:

| Building type                   |    1 |    2 |    3 |    4 |    5 |    6 |
|---------------------------------|------|------|------|------|------|------|
| `emis_domestic_compact_factors` | 0.23 | 0.28 | 0.28 | 0.26 | 0.29 | 0.29 | 
| `emis_domestic_energy_demands`  |  130 |  100 |  100 |  11  |   89 |   89 |

`emis_domestic_species_names` and `emis_domestic_species_emission_factors`: Emission factors are to be presented on a name-value pair represented by these two arrays in the `_p3d` file. Emission factors are expressed in a per terajoule (TJ) of energy consumed, and are to be presented in moles for all gas phase species, and kilograms for particulate species. In theory, the user can define up to 20 species, each represented by a name of up to 64 characters. However, only names that are found in the active chemical mechanism will have their corresponding emission factors extracted and used in the parametrization. Due to the variability species names in active chemical mechanisms, the user must provide all species names and corresponding emission factors.

The following table summarizes the recommended values for emission factors.  Speciation are required for NO<sub>x</sub> and VOC emissions for appropriation into active mechanism species:

| Heating type | CO (mol/TJ) | NO<sub>2</sub> (mol/TJ) | PM<sub>10</sub> (kg/TJ) | NO<sub>x</sub> (kg/TJ) | VOC (kg/TJ) |
|--------------------------|-------|------|-------|----|-----|
| Centralized oil          |  0.1  | 2.1  | 0.34  | 45 | 0.5 |
| Centralized gas          |  0.14 | 0.78 | 0.006 | 17 | 0.7 |
| Centralized wood pellets |  1.7  | 3.4  | 18    | 73 | 3.2 |
| Centralized woodchips    |  1.6  | 4.2  | 27    | 91 | 1.8 |
| Centralized logs         |  8.3  | 3.9  | 40    | 84 | 22  |
| Wood stoves / fireplaces | 28    | 3.9  | 48    | 84 | 29  |

## Specifications under LOD 2

Input data for the domestic heating emission mode are contained in a file called `[model]_emis_domestic` in the model `INPUT` directory, where `[model]` refers to the name of the model. The file is to be stored using netCDF APIs following the format for [specific emission mode input data in LOD 2](./EMISSIONS_LOD2_spec.md).
