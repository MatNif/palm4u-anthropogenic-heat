# Emission Module

The emission module allows pollutant sources to be specified at the available level of detail (LOD) specified by the user. There are two possible methods to introduce emissions in the `_p3d` file:

1. By setting the `emis_anthropogenic` to `.TRUE.`
2. By setting the specific emission modes to `.TRUE.`

Both options can be used simultaneously and independently in a given model.

The emission module is a component of the chemistry module, in which only species defined in the active chemical mechanism used will be used. Input emission species absent from the active mechanism will be ignored.

## Levels of Detail (LOD) and Units of Input

Depending on the model objectives or available input data, emissions can be specified at specific levels of detail (LOD). Up to three levels LODs are available:

* LOD 0 : Emissions are parameterized based on aggregated values provided by the user in the `_p3d` file
* LOD 1 : Emissions are parameterized based in part on user-defined aggregated values in the `_p3d` file, as well as temporal or spatial specific input data stored in external netCDF files
* LOD 2 : Emissions are fully specified by the user using external input data stored in external netCDF files

In legacy model specifications the three LODs are also referred to as *PARAMETERIZED*, *DEFAULT*, and  *PRE-PREOCESSED* modes respectively. This nomenclature has been deprecated.

While aggregate emissions in LOD 0 or LOD 1 are provided in units specific to each `_p3d` option, emissions in LOD 2 should be provided in SI units, that is, in kg, moles, m, and s. Further, unless otherwise indicated, reactive gas phase emission species, such as NO<sub>2</sub> or O<sub>3</sub> are to be provided in moles. On  the other hand, tracer or inert species, such as PM<sub>10</sub> or pollen (in the pollen module), are to be specified in kilograms.

## Anthropogenic Emissions Model

Emissions from the anthropogenic emissions model are specified as surface fluxes on the lowest vertical level of each domain cell column. It can be activated in the `_p3d` file by setting the option `emis_anthropogenic` to `.TRUE.`Further information on using the anthropogenic emission model can be found the in the following link:

* [Anthropogenic emissions model](./CS_model.md)

## Specific Emission Modes

Alternatively, emissions from different modes or sectors can be invoked independently into the model. This enables emission data sets from specific sectors at corresponding LODs to be prepared independently. The following emission modes are currently supported:

| Emission mode    | `_p3d` option     | Supported LOD(s) |
|------------------|-------------------|------------------|
| Biogenic         | `emis_biogenic`   | 0                |
| Domestic heating | `emis_domestic`   | 2                |
| Pollen           | `emis_pollen`     | 0                |
| Point sources    | `emis_pt_source`  | 0                |
| Traffic          | `emis_traffic`    | 2                |

Additionally, a generic mode has been implemented to accommodate emission modes that have not been implemented, or in cases where emissions from multiple modes have been assembled into a single data set. It can be activated by setting the value of the `_p3d` option `emis_generic` to `.TRUE.`.  The generic mode is only availabe in LOD 2.

As previously indicated, inputs for LODs 0 and 1 are mode-specific and will be described individually. Specification of all emission modes under LOD 2, when supported, all follow the same netCDF file format. Information pertaining to preparation of the LOD 2 data set can be found in the following link:

* [Specific emission mode input data in LOD 2](./EMISSIONS_LOD2_spec.md)

Further information for all inputs relevant to each emission mode can be found in the following links:

* [Biogenic emissions model](./BVOC_model.md)
* [Domestic heating emissions mode](./DOMESTIC_model.md)
* [Pollen model](./POLLEN_model.md)
* [Point source emissions mode](./PTSRC_model.md)
* [Traffic emissions mode](./TRAFFIC_model.md)
* [Generic emissions mode](./GENERIC_model.md)
